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
		convert << parameter.run_id << "/" << parameter.run_id << ".parameter"; 
		file_name = storage_filename_base + convert.str(); 
		parameter.SaveParameterToFile(file_name); 
	}

	/* Seed the salves; send out one unit work to each slave */
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

	// tell all the slaves to exit by sending an empty messag with 0 simulation length 
	for (int rank=1; rank<nTasks; rank++)
		MPI_Send(0, 0, MPI_INT, rank, 0, MPI_COMM_WORLD);

	parameter.number_cluster_node = nTasks;
	parameter.TraceStorageHead(storage);
        
       	convert.str(string());
        convert <<  parameter.run_id << "/" << parameter.run_id << ".parameter";
        file_name = storage_filename_base +  convert.str();
        parameter.SaveParameterToFile(file_name);

	convert.str(string()); 
	convert << parameter.run_id << "/" << parameter.run_id << ".summary";
 	file_name = storage_filename_base + convert.str(); 
        parameter.WriteSummaryFile(file_name);
	
	// Write current state
	convert.str(string());
        convert << parameter.run_id << "/" << parameter.run_id << ".current_state." << my_rank;
        file_name = storage_filename_base + convert.str();
        parameter.SaveCurrentStateToFile(file_name);
}

