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
void RunSimulation(CParameterPackage &, int, int, CModel *, CStorageHeadPthread &, const gsl_rng *); 

void slave(string storage_filename_base, CStorageHeadPthread &storage, CParameterPackage &parameter, int highest_level, bool if_resume, CModel *target, const gsl_rng* r)
{
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  
	MPI_Status status; 
	int start; 
	while (1)
	{
		MPI_Recv(&start, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
		if (status.MPI_TAG == 0)
			return; 
		else if (status.MPI_TAG < 0)
			exit(-1);

		RunSimulation(parameter, highest_level, my_rank, target, storage, r); 
	
		parameter.TraceStorageHead(storage); 
		stringstream convert; 
		convert.str(string()); 
		convert << parameter.run_id << "/" << parameter.run_id << ".current_state." << my_rank; 
		string file_name = storage_filename_base + convert.str(); 
		parameter.SaveCurrentStateToFile(file_name); 
		
		int result = 1; 
		/* Send info back */
		MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); 
	}
}

