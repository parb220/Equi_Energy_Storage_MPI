/* Test: 
 * CEquiEnergy
 * CModel
 * CMixtureModel
 * CTransitionModel
 * CSimpleGaussianModel
 * CTransitionModel_SimpleGaussian
 * CTransitionModel_Gaussian
 * CUniformModel
 * CBoundedModel
 * CEES_Head
 * CEES_Node
 * CStorageHead
 * CPutGetBin
 * CSampleIDWeight
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <mpi.h>
#include "CMixtureModel.h"
#include "CSimpleGaussianModel.h"
#include "CParameterPackage.h"
#include "CStorageHeadPthread.h"
#include "equi_energy_setup_constant.h"

using namespace std;

bool Configure_GaussianMixtureModel_File(CMixtureModel &, const string); 
void master(string, CStorageHeadPthread &, CParameterPackage &, int, bool, CModel *, const gsl_rng *); 
void slave(string, CStorageHeadPthread &, CParameterPackage &, int, bool, CModel *, const gsl_rng *); 

void usage(int arc, char **argv)
{
        cerr << "usage: " << argv[0] << endl;
        cerr << "-i <id>: id of simulation run\n";
        cerr << "-y: to continue a previous simulation run (when -i is provided)\n";
        cerr << "-d <dimension>: dimension of samples\n";
        cerr << "-f <file>: prefix of the files of the target model\n";
        cerr << "-p <probability>: probability of equi-energy jump\n";
        cerr << "-h <energy>: energy bound of the highest energy level\n";
        cerr << "-l <length>: simulation length\n";
        cerr << "-c <C factor>: c factor to determine temperature bounds according to energy bounds\n";
        cerr << "-t <number>: number of tuning times\n";
        cerr << "-b <path>: directory to store samples\n";
        cerr << "-e <level>: highest energy level to run simulation\n";
        cerr << "? this message\n";
}

int main(int argc, char ** argv)
{
	/* Initialize MPI */
	MPI_Init(&argc, &argv); 
	int my_rank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	
	/* Initialize the random_number_generator */
	const gsl_rng_type *T; 
	gsl_rng *r;
	gsl_rng_env_setup(); 
	T = gsl_rng_default; 
	r = gsl_rng_alloc(T); 
	gsl_rng_set(r, (unsigned)time(NULL)); 	
	
	/* default setting; 	*/
        string target_filename_base = string("../equi_energy_generic/gaussian_mixture_model.");
	int _run_id = time(NULL); // by default, use current time as run_id; 
        bool if_resume = false;
        string storage_filename_base = string("/home/f1hxw01/equal_energy_hw/equi_energy_storage_data/");
        int _data_dimension = DATA_DIMENSION;
        double _pee = PEE;
        double _h_k_1 = HK_1;
        int _simulation_length = SIMULATION_LENGTH;
        double _c_factor = C;
        double _mh_target_acc = MH_TARGET_ACC;
        double _energy_level_tuning_max_time = ENERGY_LEVEL_TUNING_MAX_TIME;
	double highest_level = -1; 

	/* parse command line options */
	int opt;
        while ( (opt = getopt(argc, argv, "i:yd:f:p:h:l:c:t:b:e:?")) != -1)
        {
                switch (opt)
                {
			case 'i':
                                _run_id = atoi(optarg); break;
			case 'y':
                                if_resume = true; break;
			case 'd':
                                _data_dimension = atoi(optarg); break;
                        case 'f':
				target_filename_base = string(optarg); break;
			case 'p':
                                _pee = atof(optarg); break;
			case 'h':
                                _h_k_1 = atof(optarg); break;
			case 'l':
                                _simulation_length = atoi(optarg); break;
			case 'c':
                                _c_factor = atof(optarg); break;
			case 't':
                                _energy_level_tuning_max_time = atoi(optarg); break;
			case 'b':
                                storage_filename_base = string(optarg); break;
			case 'e':
				highest_level = atoi(optarg); break; 
			case '?':
                        {
                                usage(argc, argv);
                                exit(-1);
                        }
                }
        }
	// Initialize parameters
	CParameterPackage parameter;
        stringstream convert;
	string file_name; 
        if (if_resume)
        {
                convert.str(std::string());
                convert <<  _run_id  << "/" << _run_id << ".parameter";
                file_name = storage_filename_base + convert.str();
                parameter.LoadParameterFromFile(file_name);
        }
	else 
	{
		parameter.run_id = _run_id;
                parameter.get_marker = 10000;
                parameter.put_marker = 10000;
                parameter.number_energy_level = NUMBER_ENERGY_LEVEL;
                parameter.data_dimension = _data_dimension;
                parameter.number_bins = parameter.number_energy_level * parameter.number_energy_level;
                parameter.pee = _pee;
                parameter.h0 = H0;
                parameter.hk_1 = _h_k_1;
                parameter.energy_tracking_number = ENERGY_TRACKING_NUMBER;
                parameter.t0 = T0;
                parameter.c_factor = _c_factor;
                parameter.mh_target_acc = _mh_target_acc;
                parameter.initial_sigma = INITIAL_SIGMA;
                parameter.uniform_lb = 0.0;
                parameter.uniform_ub = 1.0;
                parameter.burn_in_period = BURN_IN_PERIOD;
                parameter.multiple_try_mh = MULTIPLE_TRY_MH;
                parameter.mh_tracking_length = MH_TRACKING_LENGTH;
                parameter.mh_stepsize_tuning_max_time = MH_STEPSIZE_TUNING_MAX_TIME;
                parameter.energy_level_tracking_window_length = ENERGY_LEVEL_TRACKING_WINDOW_LENGTH;
                parameter.energy_level_tuning_max_time = _energy_level_tuning_max_time;
                parameter.deposit_frequency = DEPOSIT_FREQUENCY;

                if (MH_BLOCK)
                        parameter.number_block = parameter.data_dimension;
                else
                        parameter.number_block = 1;
                parameter.SetBlock();
	        parameter.SetEnergyBound();
                parameter.SetTemperature();
                parameter.SetMHProposalScale();
	}
	parameter.simulation_length = _simulation_length; 
	if (!if_resume || highest_level <0 || highest_level >= parameter.number_energy_level)
		highest_level = parameter.number_energy_level-1;
	
	/* Initialize the target distribution as a Gaussian mixture model; */
	// Mean Sigma and Weight are stored in files
	CMixtureModel target; 
	if (!Configure_GaussianMixtureModel_File(target, target_filename_base))
	{
		cerr << "Error in configuring gaussian mixture model as the target model.\n"; 
		exit (-1);
	}
	// Storage
	CStorageHeadPthread storage(parameter.run_id, parameter.get_marker, parameter.put_marker, parameter.number_bins,storage_filename_base, my_rank);
	if (if_resume)
                storage.restore();
        else
                storage.makedir();

	// Current state
	if (if_resume)
	{
		convert.str(string()); 
		convert << parameter.run_id << "/" << parameter.run_id << ".current_state." << my_rank; 
		file_name = storage_filename_base + convert.str(); 
		if (!parameter.LoadCurrentStateFromFile(file_name) && !parameter.LoadCurrentStateFromStorage(storage, r))
			parameter.SetCurrentState(r); 
	}
	else 
		parameter.SetCurrentState(r); 
	
	// master runs burn-in and tuning and simulation 
	// slave runs simulation
	if (my_rank == 0)
		master(storage_filename_base, storage, parameter, highest_level, if_resume, &target, r); 
	else 
		slave(storage_filename_base, storage, parameter, highest_level, if_resume, &target, r); 

	gsl_rng_free(r);	
	/* Shut down MPI */
	MPI_Finalize(); 
	return 0; 
}

bool Configure_GaussianMixtureModel_File(CMixtureModel &mixture_model, const string filename_base)
{
	/*weight */
	string filename = filename_base + "weight";
	int nComponent, nDim; 			// Number of components, dimension of variables
	bool equalComponent, equalDim;		// Whether parameters for different components are thes same; and whether parameters for different dimension are the same; 
	ifstream inputFile; 
	inputFile.open(filename.data());
	if (!inputFile)
		return false;
	inputFile >> nComponent; 
	double *weight = new double[nComponent]; 
	inputFile >> equalComponent; 
	if (!equalComponent)
	{
		for (int i=0; i<nComponent; i++)
			inputFile >> weight[i]; 
	} 
	else 
	{
		inputFile >> weight[0]; 
		for (int i=1; i<nComponent; i++)
			weight[i] = weight[0];
	}
	mixture_model.SetModelNumber(nComponent); 
	mixture_model.SetWeightParameter(weight, nComponent); 
	delete [] weight; 
	inputFile.close();

	/*sigma */
	filename = filename_base + "sigma"; 
	inputFile.open(filename.data()); 
	if (!inputFile)
		return false; 
	inputFile >> nComponent >> nDim; 
	double** sigma = new double* [nComponent];
	for (int i=0; i<nComponent; i++)
		sigma[i] = new double[nDim]; 
	inputFile >> equalComponent >> equalDim; 
	if (!equalComponent && !equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				inputFile >> sigma[i][j]; 
		}
	}
	else if (equalComponent && !equalDim)
	{
		for (int j=0; j<nDim; j++)
			inputFile >> sigma[0][j]; 
		for (int i=1; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				sigma[i][j] = sigma[0][j];
		}
	}
	else if (!equalComponent && equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			inputFile >> sigma[i][0];
			for (int j=1; j<nDim; j++)
				sigma[i][j] = sigma[i][0]; 
		} 
	}
	else
	{
		inputFile >> sigma[0][0]; 
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				sigma[i][j] = sigma[0][0];
		}
	}

	inputFile.close(); 

	/*mean */
	filename = filename_base + "mean"; 
	inputFile.open(filename.data()); 
	if (!inputFile)
		return false; 
	inputFile >> nComponent >> nDim; 
	double** mean = new double* [nComponent];
	for (int i=0; i<nComponent; i++)
		mean[i] = new double[nDim]; 
	inputFile >> equalComponent >> equalDim; 
	if (!equalComponent && !equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				inputFile >> mean[i][j]; 
		}
	}
	else if (equalComponent && !equalDim)
	{
		for (int j=0; j<nDim; j++)
			inputFile >> mean[0][j]; 
		for (int i=1; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				mean[i][j] = mean[0][j];
		}
	}
	else if (!equalComponent && equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			inputFile >> sigma[i][0];
			for (int j=1; j<nDim; j++)
				mean[i][j] = mean[i][0]; 
		} 
	}
	else
	{
		inputFile >> mean[0][0]; 
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				mean[i][j] = mean[0][0];
		}
	}

	inputFile.close(); 

	mixture_model.SetDataDimension(nDim); 
	for (int i=0; i<nComponent; i++)
	{
		mixture_model.Initialize(i, new CSimpleGaussianModel(nDim, mean[i], sigma[i])); 
		delete []mean[i]; 
		delete []sigma[i]; 
	}
	mixture_model.CalculateSetParameterNumber();
	return true;
}
