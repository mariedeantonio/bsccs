/*
 * GPUCyclicCoordinateDescent.cpp
 *
 *  Created on: July, 2010
 *      Author: msuchard
 */

#include <iostream>
#include <cmath>
#include <omp.h>

#include "GPUCyclicCoordinateDescent.h"
#include "GPU/GPUInterface.h"
#include "GPU/KernelLauncherCCD.h"
#include "GPUContainer.h"


using namespace std;

GPUCyclicCoordinateDescent::GPUCyclicCoordinateDescent(int deviceNumber, InputReader* reader)
	: CyclicCoordinateDescent(reader) {

#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tEntering GPUCylicCoordinateDescent::constructor\n");
#endif

	cout << "Running GPU version" << endl;
#ifdef MULTI_GPU

	firstGPU = new GPUInterface;
	gpuDeviceCount2 = 0;
	if (firstGPU->Initialize()) {
		gpuDeviceCount2 = firstGPU->GetDeviceCount();
		cout << "Number of GPU devices found: " << gpuDeviceCount2 << endl;
	} else {
		cerr << "Unable to initialize CUDA driver library!" << endl;
		exit(-1);
	}

	omp_set_num_threads(gpuDeviceCount2);
	partition = new GPUContainer*[gpuDeviceCount2];

	vector<int> rowOffsets(N+1);

	int offset = 0;
	int currentPid = -1;
	for (int i = 0; i < K; i++) {
		int thisPid = hPid[i];
		if (thisPid != currentPid) {
			rowOffsets[thisPid] = offset;
			currentPid = thisPid;
		}
		offset++;
	}
	rowOffsets[N] = offset;
#pragma omp parallel
	{
		int g = omp_get_thread_num();
		
		partition[g] = new GPUContainer;
		
		partition[g]->gpu = new GPUInterface;

		int numberChars = 80;
		std::vector<char> charVector2 (numberChars);
		partition[g]->gpu->GetDeviceName(g, &charVector2[0], numberChars);
		string deviceNameString(charVector2.begin(), charVector2.end());
		partition[g]->gpu->GetDeviceDescription(0, &charVector2[0]);
		string deviceDescrString(charVector2.begin(), charVector2.end());
		cout << "Using "<< deviceNameString << ": " << deviceDescrString << endl;
		partition[g]->kernels = new KernelLauncherCCD(partition[g]->gpu);
		partition[g]->gpu->SetDevice(g, partition[g]->kernels->getKernelsString());
		partition[g]->kernels->LoadKernels();
#pragma omp barrier
	}
		gpu = partition[0]->gpu;
                dBeta = gpu->AllocateRealMemory(J);
                gpu->MemcpyHostToDevice(dBeta, hBeta, sizeof(REAL) * J);

	for (int g = 0; g < gpuDeviceCount2; g++)
	{
		partition[g]->hReader = reader;
		partition[g]->drugs = J/gpuDeviceCount2;
		partition[g]->patients = N/gpuDeviceCount2;
		if (g == (gpuDeviceCount2-1))
		{
			partition[g]->drugs+=(J%gpuDeviceCount2);
			partition[g]->patients+=(N%gpuDeviceCount2);
		}
		partition[g]->prev_exposures = 0;
		partition[g]->prev_patients = 0;
		for (int i = 0; i < g; i++)
		{
			partition[g]->prev_exposures+= partition[i]->exposures;
			partition[g]->prev_patients += partition[i]->patients;
		}
		if (g == (gpuDeviceCount2-1))
			partition[g]->exposures = K-partition[g]->prev_exposures;
		else
			partition[g]->exposures = rowOffsets[partition[g]->patients + partition[g]->prev_patients]-partition[g]->prev_exposures;
		
		partition[g]->dXI = (GPUPtr*) malloc(J * sizeof(GPUPtr));
		partition[g]->columnLength = new int[J];
		for (int j = 0; j < J; j++) {
			partition[g]->columnLength[j] = hXI->getNumberOfEntries(j);
			if (partition[g]->columnLength[j] != 0) {
				int temp_length = partition[g]->columnLength[j];
				partition[g]->columnLength[j] = 0;
				int* starting_value = 0;
				
				vector<int> temp_column;
				
				for (int i = 0; i < temp_length; i++)
				{
					if ((hXI->getCompressedColumnVector(j))[i] >= partition[g]->prev_exposures
							&& (hXI->getCompressedColumnVector(j))[i] < (partition[g]->exposures + partition[g]->prev_exposures))
					{
						/*if (partition[g]->columnLength[j] == 0)
						{
							starting_value = &((hXI->getCompressedColumnVector(j))[i]);
						}*/
						temp_column.push_back((hXI->getCompressedColumnVector(j))[i] - partition[g]->prev_exposures);
						partition[g]->columnLength[j]++;
						
					}
				}
				if (partition[g]->columnLength[j]!=0)
				{
					partition[g]->dXI[j] = partition[g]->gpu->AllocateRealMemory(partition[g]->columnLength[j]);
					partition[g]->gpu->MemcpyHostToDevice(partition[g]->dXI[j], &temp_column[0],
							sizeof(int) * partition[g]->columnLength[j]);

				}
				else
				{
					partition[g]->dXI[j] = 0;
				}

			}
		}
	//	dXColumnLength = gpu->AllocateIntMemory(K);
	//	gpu->MemcpyHostToDevice(dXColumnLength, &columnLength[0], sizeof(int) * K);
		//	cerr << "Memory allocate 1" << endl;
		// Allocate GPU memory for X and beta
		

//		partition[g]->dBeta = partition[g]->gpu->AllocateRealMemory(partition[g]->drugs);
//		partition[g]->gpu->MemcpyHostToDevice(partition[g]->dBeta, (hBeta+(g*(int)(J/gpuDeviceCount2))), sizeof(REAL) * partition[g]->drugs);
		partition[g]->dXBeta = partition[g]->gpu->AllocateRealMemory(partition[g]->exposures);
		partition[g]->gpu->MemcpyHostToDevice(partition[g]->dXBeta, (hXBeta+partition[g]->prev_exposures), sizeof(REAL) * partition[g]->exposures);
		//	cerr << "Memory allocate 2" << endl;
		// Allocate GPU memory for integer vectors
		partition[g]->dOffs = partition[g]->gpu->AllocateIntMemory(partition[g]->exposures);
		partition[g]->gpu->MemcpyHostToDevice(partition[g]->dOffs, (hOffs+partition[g]->prev_exposures), sizeof(int) * partition[g]->exposures);
		partition[g]->dEta = partition[g]->gpu->AllocateIntMemory(partition[g]->exposures);
		partition[g]->gpu->MemcpyHostToDevice(partition[g]->dEta, (hEta+partition[g]->prev_exposures), sizeof(int) * partition[g]->exposures);
		partition[g]->dNEvents = partition[g]->gpu->AllocateIntMemory(partition[g]->patients);
	//	gpu->MemcpyHostToDevice(dNEvents, hNEvents, sizeof(int) * N); // Moved to computeNEvents
	//	dPid = gpu->AllocateIntMemory()
		//	cerr << "Memory allocate 3" << endl;
		// Allocate GPU memory for intermediate calculations
		partition[g]->dOffsExpXBeta = partition[g]->gpu->AllocateRealMemory(partition[g]->exposures);
		partition[g]->dXOffsExpXBeta = partition[g]->gpu->AllocateRealMemory(partition[g]->exposures);
		partition[g]->dDenomPid = partition[g]->gpu->AllocateRealMemory(partition[g]->patients);
		//	real *tmp = (real *)calloc(sizeof(real), N);
	//	gpu->MemcpyHostToDevice(dDenomPid, tmp, sizeof(real) * N);
	//	free(tmp);
	//	tmp = (real *)calloc(sizeof(real), K);
	//	gpu->MemcpyHostToDevice(dOffsExpXBeta, tmp, sizeof(real) * K);
	//	free(tmp);
		//	cerr << "Memory allocate 4" << endl;


		partition[g]->dNumerPid = partition[g]->gpu->AllocateRealMemory(partition[g]->patients);
		partition[g]->dT1 = partition[g]->gpu->AllocateRealMemory(partition[g]->patients);
#ifdef GRADIENT_HESSIAN_GPU
#ifndef GH_REDUCTION_GPU
		partition[g]->dGradient = partition[g]->gpu->AllocateRealMemory(partition[g]->patients);
		partition[g]->dHessian = partition[g]->gpu->AllocateRealMemory(partition[g]->patients);
#else
		partition[g]->dGradient = partition[g]->gpu->AllocateRealMemory(2 * partition[g]->patients);
		partition[g]->dHessian = partition[g]->dGradient + sizeof(real) * partition[g]->patients;
		partition[g]->dReducedGradientHessian = partition[g]->gpu->AllocateRealMemory(2);
#endif
		partition[g]->hGradient = (real*) malloc(sizeof(real) * partition[g]->patients);
		partition[g]->hHessian = (real*) malloc(sizeof(real) * partition[g]->patients);
#endif

	//	cerr << "Memory allocate 5" << endl;
	// Allocate computed indices for sparse matrix operations



		partition[g]->dXFullRowOffsets = partition[g]->gpu->AllocateIntMemory(partition[g]->patients+1);

		for (int i = 0; i < partition[g]->patients+1; i++)
		{
			rowOffsets[partition[g]->prev_patients+i] -= partition[g]->prev_exposures;
		}

		partition[g]->gpu->MemcpyHostToDevice(partition[g]->dXFullRowOffsets,
				&rowOffsets[partition[g]->prev_patients], sizeof(int) * (partition[g]->patients+1));
		rowOffsets[partition[g]->prev_patients+partition[g]->patients]+=partition[g]->prev_exposures;

		//	cerr << "Memory allocate 6" << endl;


		partition[g]->dXColumnRowIndicators = (GPUPtr*) malloc(J * sizeof(GPUPtr));
	//	hColumnRowLength = (int*) malloc(J * sizeof(int));
		unsigned int maxActiveWarps = 0;
		for (int j = 0; j < J; j++) {
			const int n = partition[g]->columnLength[j];
			if (n > 0) {
				if (n / WARP_SIZE > maxActiveWarps) { // TODO May be reduce for very large entries
					maxActiveWarps = n / WARP_SIZE;
				}
				vector<int> columnRowIndicators(n);
				int prev_partition_column = 0;
				for (int k = 0; k < g; k++)
				{
					prev_partition_column+= partition[k]->columnLength[j];
					
				}
				const int* indicators = hXI->getCompressedColumnVector(j);
				for (int i = 0; i < n; i++) { // Loop through non-zero entries only
					int thisPid = hPid[indicators[i+prev_partition_column]];
					columnRowIndicators[i] = thisPid - partition[g]->prev_patients;
				}
				partition[g]->dXColumnRowIndicators[j] = partition[g]->gpu->AllocateIntMemory(n);
				partition[g]->gpu->MemcpyHostToDevice(partition[g]->dXColumnRowIndicators[j],
						&columnRowIndicators[0], sizeof(int) * n);
			} else {
				partition[g]->dXColumnRowIndicators[j] = 0;
			}
		}

/*	unsigned int maxActiveWarps = 0;
        for (int j = 0; j < J; j++) {
                const int n = hXI->getNumberOfEntries(j);
                if (n > 0) {
                        if (n / WARP_SIZE > maxActiveWarps) { // TODO May be reduce for very large entries
                                maxActiveWarps = n / WARP_SIZE;
                        }
                        vector<int> columnRowIndicators(n);
                        const int* indicators = hXI->getCompressedColumnVector(j);
                        for (int i = 0; i < n; i++) { // Loop through non-zero entries only
                                int thisPid = hPid[indicators[i]];
                                columnRowIndicators[i] = thisPid;
                        }
                        dXColumnRowIndicators[j] = gpu->AllocateIntMemory(n);
                        gpu->MemcpyHostToDevice(dXColumnRowIndicators[j],
                                        &columnRowIndicators[0], sizeof(int) * n);
                } else {
                        dXColumnRowIndicators[j] = 0;
                }
        }*/

		//	cerr << "Memory allocate 7" << endl;
		maxActiveWarps++;
		partition[g]->dTmpCooRows = partition[g]->gpu->AllocateIntMemory(maxActiveWarps);
		partition[g]->dTmpCooVals = partition[g]->gpu->AllocateRealMemory(maxActiveWarps);
		cout << "MaxActiveWarps = " << maxActiveWarps << endl;
	}
#else
	gpu = new GPUInterface;
	int gpuDeviceCount = 0;
	if (gpu->Initialize()) {
		gpuDeviceCount = gpu->GetDeviceCount();
		cout << "Number of GPU devices found: " << gpuDeviceCount << endl;
	} else {
		cerr << "Unable to initialize CUDA driver library!" << endl;
		exit(-1);
	}


	if (deviceNumber < 0 || deviceNumber >= gpuDeviceCount) {
		cerr << "Unknown device number: " << deviceNumber << endl;
		exit(-1);
	}

	int numberChars = 80;
	std::vector<char> charVector(numberChars);
	gpu->GetDeviceName(deviceNumber, &charVector[0], numberChars);
	string deviceNameString(charVector.begin(), charVector.end());
	gpu->GetDeviceDescription(deviceNumber, &charVector[0]);
	string deviceDescrString(charVector.begin(), charVector.end());
	cout << "Using " << deviceNameString << ": " << deviceDescrString << endl;

	kernels = new KernelLauncherCCD(gpu);
	gpu->SetDevice(deviceNumber, kernels->getKernelsString());
	kernels->LoadKernels();
	// GPU device is already to go!

	hReader = reader;

//	cerr << "Memory allocate 0" << endl;
	// Allocate GPU memory for X
	dXI = (GPUPtr*) malloc(J * sizeof(GPUPtr));
	vector<int> columnLength(J);
	for (int j = 0; j < J; j++) {
		columnLength[j] = hXI->getNumberOfEntries(j);
		if (columnLength[j] != 0) {
			dXI[j] = gpu->AllocateRealMemory(columnLength[j]);
			gpu->MemcpyHostToDevice(dXI[j], hXI->getCompressedColumnVector(j),
				sizeof(int) * columnLength[j]);
		}
	}
//	dXColumnLength = gpu->AllocateIntMemory(K);
//	gpu->MemcpyHostToDevice(dXColumnLength, &columnLength[0], sizeof(int) * K);

//	cerr << "Memory allocate 1" << endl;
	// Allocate GPU memory for X and beta
	dBeta = gpu->AllocateRealMemory(J);
	gpu->MemcpyHostToDevice(dBeta, hBeta, sizeof(REAL) * J);
	dXBeta = gpu->AllocateRealMemory(K);
	gpu->MemcpyHostToDevice(dXBeta, hXBeta, sizeof(REAL) * K);

//	cerr << "Memory allocate 2" << endl;
	// Allocate GPU memory for integer vectors
	dOffs = gpu->AllocateIntMemory(K);
	gpu->MemcpyHostToDevice(dOffs, hOffs, sizeof(int) * K);
	dEta = gpu->AllocateIntMemory(K);
	gpu->MemcpyHostToDevice(dEta, hEta, sizeof(int) * K);
	dNEvents = gpu->AllocateIntMemory(N);
//	gpu->MemcpyHostToDevice(dNEvents, hNEvents, sizeof(int) * N); // Moved to computeNEvents
//	dPid = gpu->AllocateIntMemory()

//	cerr << "Memory allocate 3" << endl;
	// Allocate GPU memory for intermediate calculations
	dOffsExpXBeta = gpu->AllocateRealMemory(K);
	dXOffsExpXBeta = gpu->AllocateRealMemory(K);
	dDenomPid = gpu->AllocateRealMemory(N);

//	real *tmp = (real *)calloc(sizeof(real), N);
//	gpu->MemcpyHostToDevice(dDenomPid, tmp, sizeof(real) * N);
//	free(tmp);
//	tmp = (real *)calloc(sizeof(real), K);
//	gpu->MemcpyHostToDevice(dOffsExpXBeta, tmp, sizeof(real) * K);
//	free(tmp);

//	cerr << "Memory allocate 4" << endl;
	dNumerPid = gpu->AllocateRealMemory(N);
	dT1 = gpu->AllocateRealMemory(N);

#ifdef GRADIENT_HESSIAN_GPU
#ifndef GH_REDUCTION_GPU
	dGradient = gpu->AllocateRealMemory(N);
	dHessian = gpu->AllocateRealMemory(N);
#else
	dGradient = gpu->AllocateRealMemory(2 * N);
	dHessian = dGradient + sizeof(real) * N;
	dReducedGradientHessian = gpu->AllocateRealMemory(2);
#endif
	hGradient = (real*) malloc(sizeof(real) * N);
	hHessian = (real*) malloc(sizeof(real) * N);
#endif

//	cerr << "Memory allocate 5" << endl;
	// Allocate computed indices for sparse matrix operations
	dXFullRowOffsets = gpu->AllocateIntMemory(N+1);
	vector<int> rowOffsets(N + 1);
	int offset = 0;
	int currentPid = -1;
	for (int i = 0; i < K; i++) {
		int thisPid = hPid[i];
		if (thisPid != currentPid) {
			rowOffsets[thisPid] = offset;
			currentPid = thisPid;
		}
		offset++;
	}
	rowOffsets[N] = offset;
	gpu->MemcpyHostToDevice(dXFullRowOffsets, &rowOffsets[0], sizeof(int) * (N + 1));

//	cerr << "Memory allocate 6" << endl;
	dXColumnRowIndicators = (GPUPtr*) malloc(J * sizeof(GPUPtr));
//	hColumnRowLength = (int*) malloc(J * sizeof(int));
	unsigned int maxActiveWarps = 0;
	for (int j = 0; j < J; j++) {
		const int n = hXI->getNumberOfEntries(j);
		if (n > 0) {
			if (n / WARP_SIZE > maxActiveWarps) { // TODO May be reduce for very large entries
				maxActiveWarps = n / WARP_SIZE;
			}
			vector<int> columnRowIndicators(n);
			const int* indicators = hXI->getCompressedColumnVector(j);
			for (int i = 0; i < n; i++) { // Loop through non-zero entries only
				int thisPid = hPid[indicators[i]];
				columnRowIndicators[i] = thisPid;
			}
			dXColumnRowIndicators[j] = gpu->AllocateIntMemory(n);
			gpu->MemcpyHostToDevice(dXColumnRowIndicators[j],
					&columnRowIndicators[0], sizeof(int) * n);
		} else {
			dXColumnRowIndicators[j] = 0;
		}
	}

//	cerr << "Memory allocate 7" << endl;
	maxActiveWarps++;
	dTmpCooRows = gpu->AllocateIntMemory(maxActiveWarps);
	dTmpCooVals = gpu->AllocateRealMemory(maxActiveWarps);

	cout << "MaxActiveWarps = " << maxActiveWarps << endl;

//	cerr << "Done with GPU memory allocation." << endl;

	//delete hReader;

	//hReader = reader; // Keep a local copy
#endif
	computeRemainingStatistics();

#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tLeaving GPUCylicCoordinateDescent::constructor\n");
#endif
}

GPUCyclicCoordinateDescent::~GPUCyclicCoordinateDescent() {
#ifdef MULTI_GPU
	//	cerr << "1" << endl;

//#pragma omp parallel
	for (int g = 0; g < gpuDeviceCount2; g++)
	{
//		unsigned int g = omp_get_thread_num();
		for (int j = 0; j < J; j++) {
			if (partition[g]->columnLength[j] > 0) {
				partition[g]->gpu->FreeMemory(partition[g]->dXI[j]);
			}
		}
		free(partition[g]->dXI);

	//	cerr << "2" << endl;
	//	partition[g]->gpu->FreeMemory(partition[g]->dBeta);
		partition[g]->gpu->FreeMemory(partition[g]->dXBeta);

		partition[g]->gpu->FreeMemory(partition[g]->dOffs);
		partition[g]->gpu->FreeMemory(partition[g]->dEta);
		partition[g]->gpu->FreeMemory(partition[g]->dNEvents);
	//	gpu->FreeMemory(dPid);
		partition[g]->gpu->FreeMemory(partition[g]->dXFullRowOffsets);

	//	cerr << "3" << endl;
		partition[g]->gpu->FreeMemory(partition[g]->dOffsExpXBeta);
		partition[g]->gpu->FreeMemory(partition[g]->dXOffsExpXBeta);
		partition[g]->gpu->FreeMemory(partition[g]->dDenomPid);
		partition[g]->gpu->FreeMemory(partition[g]->dNumerPid);
		partition[g]->gpu->FreeMemory(partition[g]->dT1);

#ifdef GRADIENT_HESSIAN_GPU
		partition[g]->gpu->FreeMemory(partition[g]->dGradient);
#ifndef GH_REDUCTION_GPU
		partition[g]->gpu->FreeMemory(partition[g]->dHessian);
#else
		partition[g]->gpu->FreeMemory(partition[g]->dReducedGradientHessian);
#endif
#endif

	//	cerr << "4" << endl;
		for (int j = 0; j < J; j++) {
			if (partition[g]->dXColumnRowIndicators[j]) {
		//		gpu->FreeMemory(dXColumnRowIndicators[j]); // TODO Causing error under Linux
			}
		}
		free(partition[g]->dXColumnRowIndicators);

	//	cerr << "5" << endl;
		partition[g]->gpu->FreeMemory(partition[g]->dTmpCooRows);
		partition[g]->gpu->FreeMemory(partition[g]->dTmpCooVals);

	//	cerr << "6" << endl;
		delete partition[g]->gpu;
		delete partition[g];
//#pragma omp barrier
	}
	delete partition;
	delete firstGPU;

#else
//	cerr << "1" << endl;
	for (int j = 0; j < J; j++) {
		if (hXI->getNumberOfEntries(j) > 0) {
			gpu->FreeMemory(dXI[j]);
		}
	}
	free(dXI);

//	cerr << "2" << endl;
	gpu->FreeMemory(dBeta);
	gpu->FreeMemory(dXBeta);

	gpu->FreeMemory(dOffs);
	gpu->FreeMemory(dEta);
	gpu->FreeMemory(dNEvents);
//	gpu->FreeMemory(dPid);
	gpu->FreeMemory(dXFullRowOffsets);

//	cerr << "3" << endl;
	gpu->FreeMemory(dOffsExpXBeta);
	gpu->FreeMemory(dXOffsExpXBeta);
	gpu->FreeMemory(dDenomPid);
	gpu->FreeMemory(dNumerPid);
	gpu->FreeMemory(dT1);

#ifdef GRADIENT_HESSIAN_GPU
	gpu->FreeMemory(dGradient);
#ifndef GH_REDUCTION_GPU
	gpu->FreeMemory(dHessian);
#else
	gpu->FreeMemory(dReducedGradientHessian);
#endif
#endif

//	cerr << "4" << endl;
	for (int j = 0; j < J; j++) {
		if (dXColumnRowIndicators[j]) {
	//		gpu->FreeMemory(dXColumnRowIndicators[j]); // TODO Causing error under Linux
		}
	}
	free(dXColumnRowIndicators);

//	cerr << "5" << endl;
	gpu->FreeMemory(dTmpCooRows);
	gpu->FreeMemory(dTmpCooVals);

//	cerr << "6" << endl;
	delete gpu;
//	cerr << "7" << endl;
#endif
}

void GPUCyclicCoordinateDescent::resetBeta(void) {

	CyclicCoordinateDescent::resetBeta();

	gpu->MemcpyHostToDevice(dBeta, hBeta, sizeof(REAL) * J);
#ifdef MULTI_GPU
#pragma omp parallel
        {
                unsigned int g = omp_get_thread_num();
                partition[g]->gpu->MemcpyHostToDevice(partition[g]->dXBeta, hXBeta+partition[g]->prev_exposures, sizeof(REAL) * partition[g]->exposures);
#pragma omp barrier
        }
#else

	gpu->MemcpyHostToDevice(dXBeta, hXBeta, sizeof(REAL) * K);
#endif
}

void GPUCyclicCoordinateDescent::computeNEvents(void) {
	CyclicCoordinateDescent::computeNEvents();
#ifdef MULTI_GPU
#pragma omp parallel
//	for (int g = 0; g < gpuDeviceCount2; g++)
	{
		unsigned int g = omp_get_thread_num();
		partition[g]->gpu->MemcpyHostToDevice(partition[g]->dNEvents,
				hNEvents+partition[g]->prev_patients, sizeof(int) * partition[g]->patients);
#pragma omp barrier
	}
#else

	gpu->MemcpyHostToDevice(dNEvents, hNEvents, sizeof(int) * N);
#endif
}

double GPUCyclicCoordinateDescent::getObjectiveFunction(void) {
//	return getLogLikelihood() + getLogPrior();
#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tEntering GPUCylicCoordinateDescent::getObjectiveFunction\n");
#endif
#ifdef MULTI_GPU
#pragma omp parallel
//    for (int g = 0; g < gpuDeviceCount2; g++)
	{
		unsigned int g = omp_get_thread_num();
    	partition[g]->gpu->MemcpyDeviceToHost(hXBeta+partition[g]->prev_exposures,
    			partition[g]->dXBeta, sizeof(real) * partition[g]->exposures);
#pragma omp barrier
    }
#else
	gpu->MemcpyDeviceToHost(hXBeta, dXBeta, sizeof(real) * K);
#endif
	double criterion = 0;
	for (int i = 0; i < K; i++) {
		criterion += hXBeta[i] * hEta[i];
	}

#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tLeaving GPUCylicCoordinateDescent::getObjectiveFunction\n");
#endif

    return criterion;
}

double GPUCyclicCoordinateDescent::computeZhangOlesConvergenceCriterion(void) {
#ifdef MULTI_GPU
#pragma omp parallel
//	for (int g = 0; g < gpuDeviceCount2; g++)
	{
		unsigned int g = omp_get_thread_num();
		partition[g]->gpu->MemcpyDeviceToHost(hXBeta+partition[g]->prev_exposures,
				partition[g]->dXBeta, sizeof(real)*partition[g]->exposures);
#pragma omp barrier
	}
#else

	gpu->MemcpyDeviceToHost(hXBeta, dXBeta, sizeof(real) * K);
#endif
	// TODO Could do reduction on GPU
	return CyclicCoordinateDescent::computeZhangOlesConvergenceCriterion();
}

void GPUCyclicCoordinateDescent::updateXBeta(double delta, int index) {

#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tEntering GPUCylicCoordinateDescent::updateXBeta\n");
#endif
	// Separate function for benchmarking
	hBeta[index] += delta;
	REAL gpuBeta = hBeta[index];
	gpu->MemcpyHostToDevice(dBeta + sizeof(REAL) * index, &gpuBeta, sizeof(REAL));

#ifdef DP_DEBUG
	fprintf(stderr,"uXB ");
	gpu->PrintfDeviceVector(dBeta, 10);
#endif
#ifdef MULTI_GPU
#pragma omp parallel
//	for (int g = 0; g < gpuDeviceCount2; g++)
	{
		int g = omp_get_thread_num();
		const int n = partition[g]->columnLength[index];
		if (n > 0)
		{
			partition[g]->kernels->updateXBeta(partition[g]->dXBeta,
					partition[g]->dXI[index], n, delta);
			partition[g]->gpu->Synchronize();
		}

#pragma omp barrier
	}


#else

	const int n = hXI->getNumberOfEntries(index);

	// OLD
#ifdef REPLICATE_ON_CPU
	const int* indicators = hXI->getCompressedColumnVector(index);
	for (int i = 0; i < n; i++) { // Loop through non-zero entries only
		const int k = indicators[i];
		hXBeta[k] += delta;
	}
#endif
	// NEW
	kernels->updateXBeta(dXBeta, dXI[index], n, delta);

#endif
#ifdef DP_DEBUG
	fprintf(stderr,"uXB ");
	gpu->PrintfDeviceVector(dXBeta, 10);
//	exit(-1);
#endif


#ifndef MULTI_GPU
#ifdef PROFILE_GPU
	gpu->Synchronize();
#endif
#endif

#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tLeaving  GPUCylicCoordinateDescent::updateXBeta\n");
#endif
}

void GPUCyclicCoordinateDescent::computeRemainingStatistics(void) {

#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tEntering  GPUCylicCoordinateDescent::computeRemainingStatistics\n");
#endif

#ifdef REPLICATE_ON_CPU
	// OLD
	zeroVector(denomPid, N);
	for (int i = 0; i < K; i++) {
		offsExpXBeta[i] = hOffs[i] * exp(hXBeta[i]);
		denomPid[hPid[i]] += offsExpXBeta[i];
	}
#endif

	// NEW
#ifdef MULTI_GPU
#pragma omp parallel
//	for (int g = 0; g < gpuDeviceCount2; g++)
	{
		unsigned int g = omp_get_thread_num();
		partition[g]->kernels->computeIntermediates(
				partition[g]->dOffsExpXBeta, //exposures
				partition[g]->dDenomPid,  //patients
				partition[g]->dOffs,  //exposures
				partition[g]->dXBeta, //exposures
				partition[g]->dXFullRowOffsets, //patients
				partition[g]->exposures,
				partition[g]->patients);
		partition[g]->gpu->Synchronize();
//		partition[g]->gpu->MemcpyDeviceToHost(denomPid+partition[g]->prev_patients, partition[g]->dDenomPid, sizeof(real)*partition[g]->patients);
	}	
#pragma omp barrier
#else

	kernels->computeIntermediates(dOffsExpXBeta, dDenomPid, dOffs, dXBeta, dXFullRowOffsets, K, N);
#endif
#ifdef MIN_GPU
	gpu->MemcpyDeviceToHost(offsExpXBeta, dOffsExpXBeta, sizeof(real) * K);
	zeroVector(denomPid, N);
	for (int i = 0; i < K; i++) {
		denomPid[hPid[i]] += offsExpXBeta[i];
	}
#else
#ifdef REDUCE_ROW_GPU
//	gpu->MemcpyDeviceToHost(denomPid, dDenomPid, sizeof(real) * N); // TODO Move to compute likelihood
#else
	gpu->MemcpyDeviceToHost(offsExpXBeta, dOffsExpXBeta, sizeof(real) * K); // TODO Remove
	gpu->MemcpyDeviceToHost(denomPid, dDenomPid, sizeof(real) * N); // TODO Move to compute likelihood
#endif
//	int I = 100;
//	cout << "Check offsExpXBeta" << endl;
//	gpu->PrintfDeviceVector(dOffsExpXBeta, I);
//	for (int i = 0; i < I; i++) {
//		cout << " " << offsExpXBeta[i];
//	}
//	cout << endl << endl;
//
//
//	cout << "Check Denom" << endl;
//
//	gpu->PrintfDeviceVector(dDenomPid, I);
//	for (int i = 0; i < I; i++) {
//		cout << " " << denomPid[i];
//	}
//	cout << endl << endl;
//	exit(0);
#endif // MIN_GPU

//	exit(0);

	sufficientStatisticsKnown = true;
#ifndef MULTI_GPU
#ifdef PROFILE_GPU
	gpu->Synchronize();
#endif
#endif


#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tLeaving  GPUCylicCoordinateDescent::computeRemainingStatistics\n");
#endif

}

void GPUCyclicCoordinateDescent::computeRatiosForGradientAndHessian(int index) {
#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tEntering GPUCylicCoordinateDescent::computeRatiosForGradientAndHessian\n");
#endif

#ifndef REDUCE_ROW_GPU
	// OLD
	zeroVector(numerPid, N);
	const int* indicators = hXI->getCompressedColumnVector(index);
	const int n = hXI->getNumberOfEntries(index);
	for (int i = 0; i < n; i++) { // Loop through non-zero entries only
		const int k = indicators[i];
		numerPid[hPid[k]] += offsExpXBeta[k];
	}

	for (int i = 0; i < N; i++) {
		t1[i] = numerPid[i] / denomPid[i];
	}
#else
	// NEW
#ifdef MULTI_GPU

#pragma omp parallel
//	for (int g = 0; g < gpuDeviceCount2; g++)
	{
			unsigned int g = omp_get_thread_num();
			const int n = partition[g]->columnLength[index];
			partition[g]->kernels->computeDerivatives(
			partition[g]->dNumerPid,
			partition[g]->dXColumnRowIndicators[index],
			partition[g]->dXI[index],
			partition[g]->dOffsExpXBeta,
			partition[g]->dDenomPid,
			partition[g]->dT1,
#ifdef GRADIENT_HESSIAN_GPU
			partition[g]->dGradient,
			partition[g]->dHessian,
			partition[g]->dNEvents,
#endif
			n,
			partition[g]->patients,
			partition[g]->dTmpCooRows, partition[g]->dTmpCooVals);

			
			
#pragma omp barrier
	}

#else
	const int n = hXI->getNumberOfEntries(index);
			kernels->computeDerivatives(
			dNumerPid,
			dXColumnRowIndicators[index],
			dXI[index],
			dOffsExpXBeta,
			dDenomPid,
			dT1,
#ifdef GRADIENT_HESSIAN_GPU
			dGradient,
			dHessian,
			dNEvents,
#endif
			n,
			N,
			dTmpCooRows, dTmpCooVals);
#endif
#endif
#ifndef GRADIENT_HESSIAN_GPU
	gpu->MemcpyDeviceToHost(t1, dT1, sizeof(REAL) * N);
#endif


#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tLeaving GPUCylicCoordinateDescent::computeRatiosForGradientAndHessian\n");
#endif
}

void GPUCyclicCoordinateDescent::getDenominators(void) {
#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tEntering GPUCylicCoordinateDescent::getDenominators\n");
#endif
#ifdef MULTI_GPU
#pragma omp parallel
//    for (int g = 0; g < gpuDeviceCount2; g++)
    {
		unsigned int g = omp_get_thread_num();
    	partition[g]->gpu->MemcpyDeviceToHost(denomPid+(partition[g]->prev_patients),
    			partition[g]->dDenomPid, sizeof(real)*partition[g]->patients);
#pragma omp barrier
    }
#else
	gpu->MemcpyDeviceToHost(denomPid, dDenomPid, sizeof(real) * N);
#endif
#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tLeaving GPUCylicCoordinateDescent::getDenominators\n");
#endif
}

unsigned int nextPow2( unsigned int x) {
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

void GPUCyclicCoordinateDescent::computeGradientAndHession(int index, double *ogradient,
		double *ohessian) {
#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tEntering GPUCylicCoordinateDescent::computeGradientAndHessian\n");
#endif
	real gradient = 0;
	real hessian = 0;
#ifdef COKI_REDUCTION

	unsigned int threads = (N < 512) ? nextPow2((N + 1)/ 2) : 256;
	unsigned int blocks = (N + (threads * 2 - 1)) / (threads * 2);
	blocks = (64 < blocks) ? 64 : blocks;

	kernels->reduceSum(dGradient, dGradient, N, blocks, threads);
	kernels->reduceSum(dHessian, dGradient+(blocks*sizeof(real)), N, blocks, threads);
	gpu-> MemcpyDeviceToHost(hGradient, dGradient, 2*blocks*sizeof(real));

	for (int i = 0; i < blocks; i++)
	{
		gradient += hGradient[i];
		hessian += hGradient[i+blocks];
	}

#else
	real tmp[2];
#ifdef MULTI_GPU

#pragma omp parallel
//	for (int g = 0; g < gpuDeviceCount2; g++)
	{
		unsigned int g = omp_get_thread_num();
		partition[g]->kernels->reduceTwo(partition[g]->dReducedGradientHessian, partition[g]->dGradient, partition[g]->patients);
		partition[g]->gpu->MemcpyDeviceToHost(partition[g]->hGradient, partition[g]->dReducedGradientHessian, sizeof(real)*2);
//		gradient += tmp[0];
//		hessian += tmp[1];
#pragma omp barrier
	}
	for (int g = 0; g < gpuDeviceCount2; g++)
	{
		gradient+= partition[g]->hGradient[0];
		hessian+= partition[g]->hGradient[1];
	}

#else

	kernels->reduceTwo(dReducedGradientHessian, dGradient, N);
	gpu->MemcpyDeviceToHost(tmp, dReducedGradientHessian, sizeof(real) * 2);


	gradient = tmp[0];
	hessian = tmp[1];
#endif
#endif

	gradient -= hXjEta[index];
	*ogradient = (double) gradient;
	*ohessian = (double) hessian;

#ifdef DP_DEBUG
	fprintf(stderr,"%5.3f %5.3f\n", *ogradient, *ohessian);
#endif

#ifdef GPU_DEBUG_FLOW
    fprintf(stderr, "\t\t\tLeaving GPUCylicCoordinateDescent::computeGradientAndHessian\n");
#endif
}

