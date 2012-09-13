#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <pthread.h>
#include <mpi.h>
#include "equi_energy_setup_constant.h"
#include "CMixtureModel.h"
#include "CEES_Pthread.h"
#include "CStorageHeadPthread.h"
#include "CParameterPackage.h"

using namespace std;

void *simulation(void*); 

void slave(int argc, char **argv, const CModel *target, const gsl_rng* r)
{
	int _run_id, _simulation_length, highest_level;  
	string storage_filename_base;
	highest_level = -1; 
	MPI_Status status; 
	while (1)
	{
		MPI_Recv(&_run_id, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
		if (status.MPI_TAG > 0)
			_simulation_length = status.MPI_TAG; 
		else if (status.MPI_TAG == 0)
			return; 
		else 
			exit(-1); 
		
		/* parse command line options */
		int opt;
        	while ( (opt = getopt(argc, argv, "i:yd:f:p:h:l:c:t:b:e:?")) != -1)
        	{
                	switch (opt)
                	{
				case 'b':
                                	storage_filename_base = string(optarg); break;
				case 'e':
					highest_level = atoi(optarg); break; 
                	}
		}

		// Initialize parameters
		CParameterPackage parameter;
		stringstream convert; 
                convert.str(std::string());
                convert << _run_id << ".parameter";
                string file_name = storage_filename_base + convert.str();
                parameter.LoadParameterFromFile(file_name);
		parameter.simulation_length = _simulation_length; 
		if (highest_level < 0 || highest_level >= parameter.number_energy_level)
			highest_level = parameter.number_energy_level-1;

		int my_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  
	
		/* Initialize Storage  */
		CStorageHeadPthread storage(parameter.run_id, parameter.get_marker, parameter.put_marker, parameter.number_bins,storage_filename_base, my_rank); 
		storage.restore(parameter); 
		
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
	
		/* Pthread */
		pthread_t *thread = new pthread_t[parameter.number_energy_level];
		temp_buffer_float = new double [parameter.data_dimension]; 
		
		int dim_cum_sum; 
		for (int i=parameter.number_energy_level-1; i>=0; i--)
		{
			if (my_rank >= parameter.number_cluster_node)
				simulator[i].Initialize(); 
			else 
			{
				CSampleIDWeight x_start = parameter.GetCurrentState(i, my_rank);
				simulator[i].Initialize(x_start); 
			}
			parameter.GetMHProposalScale(i, temp_buffer_float, parameter.data_dimension); 
			dim_cum_sum = 0; 
			for (int iBlock=0; iBlock<parameter.number_block; iBlock++)
			{
				simulator[i].SetProposal(new CTransitionModel_SimpleGaussian(parameter.GetBlockSize(iBlock), temp_buffer_float+dim_cum_sum), iBlock);
				//simulator[i].SetProposal(new CTransitionModel_Gaussian(parameter.GetBlockSize(iBlock), temp_buffer_float+dim_cum_sum), iBlock); 
				dim_cum_sum += parameter.GetBlockSize(iBlock); 
			}
		delete [] temp_buffer_float;

		// run through simulation
		cout << "Simulation for " << parameter.simulation_length << " steps.\n"; 
		vector <bool> continue_simulation(parameter.number_energy_level, true); 
		for (int i=highest_level; i>=0; i--)
        	{
			simulator[i].simulationL = parameter.simulation_length ;
               		pthread_create(&(thread[i]), NULL, simulation, (void*)(simulator+i));
        	}
        	for (int i=highest_level; i>=0; i--)
               		pthread_join(thread[i], NULL);
	
		storage.finalize(); 		// save to hard-disk of those unsaved data

		stringstream buffer(stringstream::binary|stringstream::out); 
		buffer.str(string()); 
		CSampleIDWeight x; 
		for (int i=0; i<parameter.number_energy_level; i++)
		{
			parameter.TraceSimulator(simulator[i], my_rank); 
			x = parameter.GetCurrentState(i,my_rank); 
			write(buffer, &x); 
		}
	
		/* Release dynamically allocated space */
		delete [] thread; 
		delete [] simulator; 
		
		/* Send info back */
		MPI_Send(buffer.str().c_str(), x.GetSize_Data()*parameter.number_energy_level, MPI_BYTE, 0, 0, MPI_COMM_WORLD); 
	}
}

