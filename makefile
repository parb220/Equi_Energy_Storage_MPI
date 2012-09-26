MPICXX = mpic++ 
CFLAGS := $(CFLAGS) -g -Wall 
LIBS := $(LIBS) -lpthread -lstdc++
LIBS_DIR := $(LIBS_DIR) -L/usr/lib64

EQUAL_ENERGY_HOME = /home/f1hxw01/equal_energy_hw
INCLUDE_DIR := $(INCLUDE_DIR) -I$(EQUAL_ENERGY_HOME)/include
LIBS := $(LIBS) -lgsl -lgslcblas -lm

DISTR_MODEL_DIR = $(EQUAL_ENERGY_HOME)/equi_energy_generic
DISTR_MODEL_OBJS = $(DISTR_MODEL_DIR)/CMixtureModel.o $(DISTR_MODEL_DIR)/CModel.o $(DISTR_MODEL_DIR)/CSimpleGaussianModel.o $(DISTR_MODEL_DIR)/CTransitionModel_SimpleGaussian.o $(DISTR_MODEL_DIR)/CUniformModel.o $(DISTR_MODEL_DIR)/CBoundedModel.o $(DISTR_MODEL_DIR)/AddScaledLogs.o 

SINGLE_CORE_VERSION_DIR = $(EQUAL_ENERGY_HOME)/equi_energy_storage
SINGLE_CORE_VERSION_OBJS = $(SINGLE_CORE_VERSION_DIR)/CEES_Node.o $(SINGLE_CORE_VERSION_DIR)/CPutGetBin.o $(SINGLE_CORE_VERSION_DIR)/CSampleIDWeight.o $(SINGLE_CORE_VERSION_DIR)/CStorageHead.o $(SINGLE_CORE_VERSION_DIR)/MHAdaptive.o $(SINGLE_CORE_VERSION_DIR)/CParameterPackage.o

MULTI_CORE_VERSION_DIR = $(EQUAL_ENERGY_HOME)/equi_energy_storage_pthread
MULTI_CORE_VERSION_OBJS = $(MULTI_CORE_VERSION_DIR)/CEES_Pthread.o $(MULTI_CORE_VERSION_DIR)/CStorageHeadPthread.o $(MULTI_CORE_VERSION_DIR)/pthread_initialize_simulation.o $(MULTI_CORE_VERSION_DIR)/TuneEnergyLevlesUpdateStorage.o

all:  test_GMM_mpi

test_GMM_mpi_obj = test_GMM_mpi.o mpi_master.o mpi_slave.o TuneScale_BurnIn.o RunSimulation.o $(MULTI_CORE_VERSION_OBJS) $(SINGLE_CORE_VERSION_OBJS) $(DISTR_MODEL_OBJS)

test_GMM_mpi: $(test_GMM_mpi_obj)
	$(MPICXX) $(CFLAGS) $(LIBS_DIR) $(LIBS) $(test_GMM_mpi_obj) -o test_GMM_mpi 

test_GMM_mpi.o: test_GMM_mpi.cpp
	$(MPICXX) $(CFLAGS) $(INCLUDE_DIR) -c test_GMM_mpi.cpp 

mpi_master.o: mpi_master.cpp
	$(MPICXX) $(CFLAGS) $(INCLUDE_DIR) -c mpi_master.cpp

mpi_slave.o: mpi_slave.cpp
	$(MPICXX) $(CFLAGS) $(INCLUDE_DIR) -c mpi_slave.cpp

TuneScale_BurnIn.o: TuneScale_BurnIn.cpp
	$(MPICXX) $(CFLAGS) $(INCLUDE_DIR) -c TuneScale_BurnIn.cpp

RunSimulation.o: RunSimulation.cpp
	$(MPICXX) $(CFLAGS) $(INCLUDE_DIR) -c RunSimulation.cpp

clean: 
	rm -f *.o  test_GMM_mpi 
