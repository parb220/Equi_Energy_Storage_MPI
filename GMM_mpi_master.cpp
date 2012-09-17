#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <pthread.h>
#include <mpi.h>
#include "CModel.h"
#include "CEES_Pthread.h"
#include "CStorageHeadPthread.h"
#include "CParameterPackage.h"
#include "CSampleIDWeight.h"

using namespace std;

void *initialize_simulate(void*);
void *tuning_simulation(void *); 
bool TuneEnergyLevels_UpdateStorage(CEES_Pthread *, CParameterPackage &); 
void TuneScale_BurnIn(CParameterPackage &, CModel *, CStorageHeadPthread &, const gsl_rng *); 
void RunSimulation(CParameterPackage &, int, int, CModel*, CStorageHeadPthread &, const gsl_rng *); 

void master(string storage_filename_base, CStorageHeadPthread &storage, CParameterPackage &parameter, int highest_level, bool if_resume, CModel *target, const gsl_rng *r) 
{	
	/* Find out how many processes there are in the default communicator */
	int nTasks, my_rank;   
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks); 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

	stringstream convert; 
	string file_name; 
	/* Initialize Storage  */
	if(!if_resume)
	{
		TuneScale_BurnIn(parameter, target, storage, r); 
		
		parameter.number_cluster_node = 1; 
		parameter.TraceStorageHead(storage);
		
		// Write parameters 
		convert.str(string());
		convert << parameter.run_id << ".parameter"; 
		file_name = storage_filename_base + convert.str(); 
		parameter.SaveParameterToFile(file_name); 
	}

	/* Seed the salves; send out one unit work to each slave*/
	int start; 
	for (int rank=1; rank<nTasks; rank++)
		MPI_Send(&start, 1, MPI_INT, rank, 1 , MPI_COMM_WORLD); 
	
	// Run task on master
	RunSimulation(parameter, highest_level, 0, target, storage, r); 
	
	/* receive all results from slaves */
	MPI_Status status;
	int result; 
	for (int rank=1; rank<nTasks; rank++)
		MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
	/* tell all the slaves to exit by sending an empty messag with 0 simulation length */
	for (int rank=1; rank<nTasks; rank++)
		MPI_Send(0, 0, MPI_INT, rank, 0, MPI_COMM_WORLD); 

	parameter.number_cluster_node = nTasks;
	parameter.TraceStorageHead(storage);
        
       	convert.str(string());
        convert <<  parameter.run_id << ".parameter";
        file_name = storage_filename_base +  convert.str();
        parameter.SaveParameterToFile(file_name);

	convert.str(string()); 
	convert << parameter.run_id << ".summary";
 	file_name = storage_filename_base + convert.str(); 
        parameter.WriteSummaryFile(file_name);
	
	// Write current state
	convert.str(string());
        convert << parameter.run_id << ".current_state." << my_rank;
        file_name = storage_filename_base + convert.str();
        parameter.SaveCurrentStateToFile(file_name);
}


void TuneScale_BurnIn(CParameterPackage &parameter, CModel *target, CStorageHeadPthread &storage, const gsl_rng *r)
{
	/* Initializing CEES_Pthread */ 
	CEES_Pthread::SetEnergyLevelNumber(parameter.number_energy_level); // Number of energy levels; 
	CEES_Pthread::SetEquiEnergyJumpProb(parameter.pee);		// Probability for equal energy jump
	CEES_Pthread::SetDataDimension(parameter.data_dimension); 	// Data dimension for simulation
	CEES_Pthread::ultimate_target = target;	

	double *temp_buffer_float=new double[parameter.number_energy_level]; 
	// Energy bound 
	parameter.GetEnergyBound(temp_buffer_float, parameter.number_energy_level); 
	CEES_Pthread::SetEnergyLevels(temp_buffer_float, parameter.number_energy_level);
	
	// Temperature
	parameter.GetTemperature(temp_buffer_float, parameter.number_energy_level); 
	CEES_Pthread::SetTemperatures(temp_buffer_float, parameter.number_energy_level); 
	delete [] temp_buffer_float; 

	CEES_Pthread::SetPthreadParameters(parameter.number_energy_level);	// Pthread condition and mutex
	// Block
	int *temp_buffer_int = new int[parameter.number_block]; 
	parameter.GetBlockSize(temp_buffer_int, parameter.number_block); 
	CEES_Pthread::SetBlockSize(temp_buffer_int, parameter.number_block); 	
	delete [] temp_buffer_int; 

	CEES_Pthread::storage = &storage; 

	/*  Generate K CEES_Pthread objects */
	CEES_Pthread *simulator = new CEES_Pthread[parameter.number_energy_level]; 
	for (int i=0; i<parameter.number_energy_level; i++)
	{
		simulator[i].SetID_LocalTarget(i);
		simulator[i].r = r; 	// random number generator  
		if (i < parameter.number_energy_level-1)
			simulator[i].SetHigherNodePointer(simulator+i+1);
		else 
			simulator[i].SetHigherNodePointer(NULL);
		simulator[i].mMH = parameter.multiple_try_mh; 
		simulator[i].depositFreq = parameter.deposit_frequency; 
	}

	CEES_Pthread::InitializeMinMaxEnergy(parameter.energy_tracking_number);	// For tuning energy levels based on newly identified min_energy 
	CEES_Pthread::SetTargetAcceptanceRate(parameter.mh_target_acc); 
	
	/* Pthread */
	pthread_t *thread = new pthread_t[parameter.number_energy_level];
	temp_buffer_float = new double [parameter.data_dimension]; 
	cout << "Initialize, burn in, tune/estimate MH stepsize and simulate for " << parameter.energy_level_tracking_window_length << " steps.\n"; 
	for (int i=parameter.number_energy_level-1; i>=0; i--)
	{
		simulator[i].burnInL = parameter.burn_in_period; 
		simulator[i].MHMaxTime = parameter.mh_stepsize_tuning_max_time; 
		simulator[i].MHInitialL = parameter.mh_tracking_length; 
		simulator[i].MHProposalScale = parameter.GetMHProposalScale(i); 
		simulator[i].simulationL = parameter.energy_level_tracking_window_length; 

		pthread_create(&(thread[i]), NULL, initialize_simulate, (void*)(simulator+i));
	}
	
	for (int i=parameter.number_energy_level-1; i>=0; i--)
		pthread_join(thread[i], NULL);

	int nEnergyLevelTuning = 0;
	while (nEnergyLevelTuning < parameter.energy_level_tuning_max_time)
        {       
		cout << "Energy level tuning: " << nEnergyLevelTuning << " for " << parameter.energy_level_tracking_window_length << " steps.\n"; 
		TuneEnergyLevels_UpdateStorage(simulator, parameter);
		for (int i=parameter.number_energy_level-1; i>=0; i--)
		{
			simulator[i].simulationL = parameter.energy_level_tracking_window_length; 
                       	pthread_create(&(thread[i]), NULL, tuning_simulation, (void*)(simulator+i));       
               	}
               	for (int i=parameter.number_energy_level-1; i>=0; i--)			
			pthread_join(thread[i], NULL);
		nEnergyLevelTuning ++;
	}
	
        for (int i=0; i<parameter.number_energy_level; i++)
                parameter.TraceSimulator(simulator[i]);
	delete [] simulator; 
	delete [] thread; 
}
	

