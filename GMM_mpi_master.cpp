#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <pthread.h>
#include <mpi.h>
#include "equi_energy_setup_constant.h"
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

void master(int argc, char **argv, CModel *target, const gsl_rng *r) 
{	
	/* Find out how many processes there are in the default communicator */
	int nTasks;  
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks); 
	
	/* default setting taking from equi_energy_setup_constant.h */
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
			case 'b':
                                storage_filename_base = string(optarg); break;
                        case 'd':
                                _data_dimension = atoi(optarg); break;
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
			case '?':
			{
				usage(argc, argv);
				/* tell all the slaves to exit */ 
        			for (int rank=1; rank<nTasks; rank++)
                			MPI_Send(0, 0, MPI_INT, rank, -1, MPI_COMM_WORLD);
				exit(-1); 
			}
                }
        }

	// Initialize parameters
	CParameterPackage parameter;
        if (if_resume)
        {
		/* When an old session is resumed, then parameter.number_cluster_node read from the parameter file is the number of clusters used in the last simulation. number_cluster_node of this simulation run is nTask */
		stringstream convert; 
                convert.str(std::string());
                convert << "/" << _run_id << "/" << _run_id << ".parameter";
                string file_name = storage_filename_base + convert.str();
                parameter.LoadParameterFromFile(file_name);
		
		convert.str(std::string()); 
		convert << "/" << _run_id << "/" <<_run_id << ".current_state.0";
		file_name = storage_filename_base + convert.str(); 
		parameter.LoadCurrentStateFromFile(file_name); 
        }
	else 
	{
		/* When a new session starts, only the master node does burning-in, scale tuning*/
		parameter.run_id = _run_id;
                parameter.get_marker = 50000;
                parameter.put_marker = 50000;
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
                parameter.SetCurrentState(r);
	}
	parameter.simulation_length = _simulation_length; 	

	/* Initialize Storage  */
	CStorageHeadPthread storage(parameter.run_id, parameter.get_marker, parameter.put_marker, parameter.number_bins,storage_filename_base, 0); 
	if(if_resume)
		storage.restore(); 
	else 
		storage.makedir(); 
		
	if (!if_resume)
	{
		// If this is a new simulation, tuning and burning-in is done at master node
		parameter.number_cluster_node = 1; 
		TuneScale_BurnIn(parameter, target, storage, r); 
		// parameters are saved
		stringstream convert; 
		convert.str(string());
		convert << "/" << parameter.run_id << "/" << parameter.run_id << ".parameter"; 
		string file_name = storage_filename_base + convert.str(); 
		parameter.SaveParameterToFile(file_name); 

		convert.str(string()); 
		convert << "/" << parameter.run_id << "/" << parameter.run_id << ".current_state.0"; 
		file_name = storage_filename_base + convert.str(); 
		parameter.SaveCurrentStateToFile(file_name); 
	}

	/* Seed the salves; send out one unit work to each slave*/
	for (int rank=1; rank<nTasks; rank++)
		MPI_Send(&parameter.run_id, 1, MPI_INT, rank, parameter.simulation_length, MPI_COMM_WORLD); 
	
	/* receive all results from slaves */
	// x_current for other nodes are obtained by message passing
	MPI_Status status;
	int result; 
	for (int rank=1; rank<nTasks; rank++)
		MPI_Recv(&result, 1, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 

	/* tell all the slaves to exit by sending an empty messag with 0 simulation length */
	for (int rank=1; rank<nTasks; rank++)
		MPI_Send(0, 0, MPI_INT, rank, 0, MPI_COMM_WORLD); 

	// save parameters into file
	parameter.TraceStorageHead(storage);
	parameter.number_cluster_node = nTasks;
        
	stringstream convert;
       	convert.str(string());
        convert << "/" << parameter.run_id << "/" << parameter.run_id << ".parameter";
        string file_name = storage_filename_base +  convert.str();
        parameter.SaveParameterToFile(file_name);

	convert.str(std::string()); 
	convert << parameter.run_id << ".summary";
 	file_name = storage_filename_base + convert.str(); 
        parameter.WriteSummaryFile(file_name);
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
	cout << "Initialize, burn in, tune/estimate MH stepsize and simulate for " << ENERGY_LEVEL_TRACKING_WINDOW_LENGTH << " steps.\n"; 
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
	storage.finalize(); 
	
	parameter.TraceStorageHead(storage);
        for (int i=0; i<parameter.number_energy_level; i++)
                parameter.TraceSimulator(simulator[i]);
	delete [] simulator; 
	delete [] thread; 
}
	

