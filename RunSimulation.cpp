#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <pthread.h>
#include <mpi.h>
#include "CTransitionModel_SimpleGaussian.h"
#include "CMixtureModel.h"
#include "CEES_Pthread.h"
#include "CStorageHeadPthread.h"
#include "CParameterPackage.h"

using namespace std;

void *simulation(void*);

void RunSimulation(CParameterPackage &parameter, int highest_level, int my_rank, CModel **target, CStorageHeadPthread &storage, const gsl_rng *r)
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
		simulator[i].ultimate_target = target[i]; 
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
	// Current state and scale of proposal
	for (int i=parameter.number_energy_level-1; i>=0; i--)
	{
		simulator[i].Initialize(parameter.GetCurrentState(i)); 
		
		parameter.GetMHProposalScale(i, temp_buffer_float, parameter.data_dimension); 
		dim_cum_sum = 0; 
		for (int iBlock=0; iBlock<parameter.number_block; iBlock++)
		{
			simulator[i].SetProposal(new CTransitionModel_SimpleGaussian(parameter.GetBlockSize(iBlock), temp_buffer_float+dim_cum_sum), iBlock);
			dim_cum_sum += parameter.GetBlockSize(iBlock); 
		}
	}
	delete [] temp_buffer_float;

	// run through simulation
	cout << "Simulation for " << parameter.simulation_length << " steps.\n"; 
	for (int i=highest_level; i>=0; i--)
        {
		simulator[i].simulationL = parameter.simulation_length ;
               	pthread_create(&(thread[i]), NULL, simulation, (void*)(simulator+i));
        }
        for (int i=highest_level; i>=0; i--)
              	pthread_join(thread[i], NULL);
	
	storage.finalize(); 		// save to hard-disk of those unsaved data

	for (int i=0; i<parameter.number_energy_level; i++)
		parameter.TraceSimulator(simulator[i]); 
		
	/* Release dynamically allocated space */
	delete [] thread; 
	delete [] simulator; 
}

