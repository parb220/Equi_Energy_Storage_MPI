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
void RunSimulation(CParameterPackage &, int, int, CModel*, CStorageHeadPthread &, const gsl_rng *); 

void TuneScale_BurnIn(CParameterPackage &parameter, CModel *target, CStorageHeadPthread &storage, const gsl_rng *r)
{
	/* Initializing CEES_Pthread */ 
	CEES_Pthread::SetEnergyLevelNumber(parameter.number_energy_level); // Number of energy levels; 
	CEES_Pthread::SetEquiEnergyJumpProb(parameter.pee);		// Probability for equal energy jump
	CEES_Pthread::SetDataDimension(parameter.data_dimension); 	// Data dimension for simulation
	// CEES_Pthread::ultimate_target = target;	

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
		simulator[i].ultimate_target = &(target[i]);
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
	

