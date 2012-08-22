/*
 * CyclicCoordinateDescent.cpp
 *
 *  Created on: May-June, 2010
 *      Author: msuchard
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <map>
#include <time.h>

#include "CyclicCoordinateDescent.h"
#include "InputReader.h"

#define PI	3.14159265358979323851280895940618620443274267017841339111328125

using namespace std;

namespace bsccs {
void compareIntVector(int* vec0, int* vec1, int dim, const char* name) {
	for (int i = 0; i < dim; i++) {
		if (vec0[i] != vec1[i]) {
			cerr << "Error at " << name << "[" << i << "]: ";
			cerr << vec0[i] << " != " << vec1[i] << endl;
			exit(0);
		}
	}
}

//CyclicCoordinateDescent::CyclicCoordinateDescent(
//		const char* fileNameX,
//		const char* fileNameEta,
//		const char* fileNameOffs,
//		const char* fileNameNEvents,
//		const char* fileNamePid
//	) {
//
//	hXI = new CompressedIndicatorMatrix(fileNameX);
//
//	K = hXI->getNumberOfRows();
//	J = hXI->getNumberOfColumns();
//
//	conditionId = "NA";
//
//	int lOffs;
//    hOffs = readVector<int>(fileNameOffs, &lOffs);
//
//    int lEta;
//    hEta = readVector<int>(fileNameEta, &lEta);
//
//    int lNEvents;
//    hNEvents = readVector<int>(fileNameNEvents, &lNEvents);
//
//    int lPid;
//    hPid = readVector<int>(fileNamePid, &lPid);
//
//    testDimension(lOffs, K, "hOffs");
//    testDimension(lEta, K, "hEta");
//    testDimension(lPid, K, "hPid");
//
//    N = lNEvents;
//
//    hasLog = false;
//
//    init();
//}

CyclicCoordinateDescent::CyclicCoordinateDescent() {
	// Do nothing
}

CyclicCoordinateDescent::CyclicCoordinateDescent(
			InputReader* reader
		) {
	N = reader->getNumberOfPatients();
	K = reader->getNumberOfRows();
	J = reader->getNumberOfColumns();

	classHierarchyVariance = reader->classHierarchyVariance;
	sigma2Beta = reader->sigma2Beta;


	hXI = reader;

	//cout << "CCD constructor hXI rows = " << hXI->getNumberOfRows() << endl;
	hEta = reader->getEtaVector();
	hOffs = reader->getOffsetVector();
	hNEvents = NULL;
	hPid = reader->getPidVector();

	conditionId = reader->getConditionId();

	init();
}

//CyclicCoordinateDescent::CyclicCoordinateDescent(
//		int inN,
//		CompressedIndicatorMatrix* inX,
//		int* inEta,
//		int* inOffs,
//		int* inNEvents,
//		int* inPid) :
//	N(inN), hXI(inX), hEta(inEta), hOffs(inOffs), hNEvents(inNEvents), hPid(inPid) {
//
//	K = hXI->getNumberOfRows();
//	J = hXI->getNumberOfColumns();
//
//	init();
//}

CyclicCoordinateDescent::~CyclicCoordinateDescent(void) {

	free(hPid);
	free(hNEvents);
	free(hEta);
	free(hOffs);
	
	free(hXBeta);
	free(hXBetaSave);
	free(hDelta);
	
#ifdef TEST_ROW_INDEX
	for (int j = 0; j < J; ++j) {
		if (hXColumnRowIndicators[j]) {
			free(hXColumnRowIndicators[j]);
		}
	}
	free(hXColumnRowIndicators);
#endif

	free(hXjEta);
	free(offsExpXBeta);
	free(xOffsExpXBeta);
//	free(denomPid);  // Nested in denomPid allocation
	free(numerPid);
	free(t1);
	
	if (hWeights) {
		free(hWeights);
	}
}

string CyclicCoordinateDescent::getPriorInfo() {
	stringstream priorInfo;
	if (priorType == LAPLACE) {
		priorInfo << "Laplace(";
		priorInfo << lambda;
	} else if (priorType == NORMAL) {
		priorInfo << "Normal(";
		priorInfo << sigma2Beta;
	}
	priorInfo << ")";
	return priorInfo.str();
}

void CyclicCoordinateDescent::resetBounds() {
	for (int j = 0; j < J; j++) {
		hDelta[j] = 2.0;
	}
}

void CyclicCoordinateDescent::init() {
	
	// Set parameters and statistics space
	hDelta = (bsccs::real*) malloc(J * sizeof(bsccs::real));
//	for (int j = 0; j < J; j++) {
//		hDelta[j] = 2.0;
//	}

	hBeta = (bsccs::real*) calloc(J, sizeof(bsccs::real)); // Fixed starting state
	hXBeta = (bsccs::real*) calloc(K, sizeof(bsccs::real));
	hXBetaSave = (bsccs::real*) calloc(K, sizeof(bsccs::real));

	// Set prior
	priorType = LAPLACE;
	//sigma2Beta = 1000;  tshaddox commented out, this is now an input parameter
	lambda = sqrt(2.0/sigma2Beta);
	
	// Recode patient ids
	int currentNewId = 0;
	int currentOldId = hPid[0];
	
	for(int i = 0; i < K; i++) {
		if (hPid[i] != currentOldId) {			
			currentOldId = hPid[i];
			currentNewId++;
		}
		hPid[i] = currentNewId;
	}
		
	// Init temporary variables
	offsExpXBeta = (bsccs::real*) malloc(sizeof(bsccs::real) * K);
	xOffsExpXBeta = (bsccs::real*) malloc(sizeof(bsccs::real) * K);

	// Put numer and denom in single memory block, with first entries on 16-word boundary
	int alignedLength = getAlignedLength(N);

	numerPid = (bsccs::real*) malloc(sizeof(bsccs::real) * 2 * alignedLength);
//	denomPid = (bsccs::real*) malloc(sizeof(bsccs::real) * N);
	denomPid = numerPid + alignedLength; // Nested in denomPid allocation
	t1 = (bsccs::real*) malloc(sizeof(bsccs::real) * N);
	hNEvents = (int*) malloc(sizeof(int) * N);
	hXjEta = (bsccs::real*) malloc(sizeof(bsccs::real) * J);
	hWeights = NULL;

	// Initialize XColumnRowIndicators for fast spMV

#ifdef TEST_ROW_INDEX
	hXColumnRowIndicators = (int**) malloc(J * sizeof(int*));
	unsigned int maxActiveWarps = 0;
	for (int j = 0; j < J; j++) {
		const int n = hXI->getNumberOfEntries(j);
		if (n > 0) {
			int* columnRowIndicators = (int*) malloc(n * sizeof(int));
			const int* indicators = hXI->getCompressedColumnVector(j);
			for (int i = 0; i < n; i++) { // Loop through non-zero entries only
				int thisPid = hPid[indicators[i]];
				columnRowIndicators[i] = thisPid;
			}
			hXColumnRowIndicators[j] = columnRowIndicators;
		} else {
			hXColumnRowIndicators[j] = 0;
		}
	}
#endif

	useCrossValidation = false;
	validWeights = false;
	sufficientStatisticsKnown = false;

#ifdef DEBUG	
	cerr << "Number of patients = " << N << endl;
	cerr << "Number of exposure levels = " << K << endl;
	cerr << "Number of drugs = " << J << endl;	
#endif          
}

int CyclicCoordinateDescent::getAlignedLength(int N) {
	return (N / 16) * 16 + (N % 16 == 0 ? 0 : 16);
}

void CyclicCoordinateDescent::computeNEvents() {
	zeroVector(hNEvents, N);
	if (useCrossValidation) {
		for (int i = 0; i < K; i++) {
			hNEvents[hPid[i]] += hEta[i] * int(hWeights[i]); // TODO Consider using only integer weights
		}
	} else {
		for (int i = 0; i < K; i++) {
			hNEvents[hPid[i]] += hEta[i];
		}
	}
	validWeights = true;
}

void CyclicCoordinateDescent::resetBeta(void) {
	for (int j = 0; j < J; j++) {
		hBeta[j] = 0.0;
	}
	for (int k = 0; k < K; k++) {
		hXBeta[k] = 0.0;
	}
	sufficientStatisticsKnown = false;
}

void CyclicCoordinateDescent::logResults(const char* fileName) {

	ofstream outLog(fileName);
	if (!outLog) {
		cerr << "Unable to open log file: " << fileName << endl;
		exit(-1);
	}

	InputReader* reader = dynamic_cast<InputReader*>(hXI);
	map<int, DrugIdType> drugMap = reader->getDrugNameMap();

	string sep(","); // TODO Make option

	for (int i = 0; i < J; i++) {
		outLog << conditionId << sep <<
//		i << sep <<
		drugMap[i] << sep << hBeta[i] << endl;
	}
	outLog.close();
}

double CyclicCoordinateDescent::getPredictiveLogLikelihood(bsccs::real* weights) {

	if (!sufficientStatisticsKnown) {
		computeRemainingStatistics(true);
	}

	getDenominators();

	double logLikelihood = 0;

	for (int i = 0; i < K; i++) {
		logLikelihood += hEta[i] * weights[i] * (hXBeta[i] - log(denomPid[hPid[i]]));
	}

	return logLikelihood;
}

bsccs::real CyclicCoordinateDescent::getBeta(int i) {
	if (!sufficientStatisticsKnown) {
		computeRemainingStatistics(true);
	}
	return hBeta[i];
}

double CyclicCoordinateDescent::getLogLikelihood(void) {

	if (!sufficientStatisticsKnown) {
		computeRemainingStatistics(true);
	}

	getDenominators();

	double logLikelihood = 0;

	if (useCrossValidation) {
		for (int i = 0; i < K; i++) {
			logLikelihood += hEta[i] * hXBeta[i] * hWeights[i];
		}
	} else {
		for (int i = 0; i < K; i++) {
			logLikelihood += hEta[i] * hXBeta[i];
		}
	}

	for (int i = 0; i < N; i++) {
		logLikelihood -= hNEvents[i] * log(denomPid[i]);
	}

	return logLikelihood;
}

void CyclicCoordinateDescent::getDenominators() {
	// Do nothing
}

double convertVarianceToHyperparameter(double value) {
	return sqrt(2.0 / value);
}

void CyclicCoordinateDescent::setHyperprior(double value) {
	sigma2Beta = value;
	lambda = convertVarianceToHyperparameter(value);
}

void CyclicCoordinateDescent::setPriorType(int iPriorType) {
	if (iPriorType != LAPLACE && iPriorType != NORMAL) {
		cerr << "Unknown prior type" << endl;
		exit(-1);
	}
	priorType = iPriorType;
}

void CyclicCoordinateDescent::setWeights(bsccs::real* iWeights) {

	if (iWeights == NULL) {
		std::cerr << "Turning off weights!" << std::endl;
		// Turn off weights
		useCrossValidation = false;
		validWeights = false;
		sufficientStatisticsKnown = false;
		return;
	}

	if (hWeights == NULL) {
		hWeights = (bsccs::real*) malloc(sizeof(bsccs::real) * K);
	}
	for (int i = 0; i < K; ++i) {
		hWeights[i] = iWeights[i];
	}
	useCrossValidation = true;
	validWeights = false;
	sufficientStatisticsKnown = false;
}
	
double CyclicCoordinateDescent::getLogPrior(void) {
	if (priorType == LAPLACE) {
		return J * log(0.5 * lambda) - lambda * oneNorm(hBeta, J);
	} else {
		return -0.5 * J * log(2.0 * PI * sigma2Beta) - 0.5 * twoNormSquared(hBeta, J) / sigma2Beta;		
	}
}

double CyclicCoordinateDescent::getObjectiveFunction(void) {	
//	return getLogLikelihood() + getLogPrior(); // This is LANGE
	double criterion = 0;
	if (useCrossValidation) {
		for (int i = 0; i < K; i++) {
			criterion += hXBeta[i] * hEta[i] * hWeights[i];
		}
	} else {
		for (int i = 0; i < K; i++) {
			criterion += hXBeta[i] * hEta[i];
		}
	}
	return criterion;
}

double CyclicCoordinateDescent::computeZhangOlesConvergenceCriterion(void) {
	double sumAbsDiffs = 0;
	double sumAbsResiduals = 0;
	if (useCrossValidation) {
		for (int i = 0; i < K; i++) {
			sumAbsDiffs += abs(hXBeta[i] - hXBetaSave[i]) * hEta[i] * hWeights[i];
			sumAbsResiduals += abs(hXBeta[i]) * hEta[i] * hWeights[i];
		}
	} else {
		for (int i = 0; i < K; i++) {
			sumAbsDiffs += abs(hXBeta[i] - hXBetaSave[i]) * hEta[i];
			sumAbsResiduals += abs(hXBeta[i]) * hEta[i];
		}
	}
	return sumAbsDiffs / (1.0 + sumAbsResiduals);
}

void CyclicCoordinateDescent::saveXBeta(void) {
	memcpy(hXBetaSave, hXBeta, K * sizeof(bsccs::real));
}

void CyclicCoordinateDescent::update(
		int maxIterations,
		int convergenceType,
		double epsilon
		) {

	if (convergenceType != LANGE && convergenceType != ZHANG_OLES) {
		cerr << "Unknown convergence criterion" << endl;
		exit(-1);
	}

	if (!validWeights) {
		computeXjEta();
		computeNEvents();
	}

	if (!sufficientStatisticsKnown) {
		computeRemainingStatistics(true);
	}

	resetBounds();

	bool done = false;
	int iteration = 0;
	double lastObjFunc;

	if (convergenceType == LANGE) {
		lastObjFunc = getObjectiveFunction();
	} else { // ZHANG_OLES
		saveXBeta();
	}
	
	while (!done) {
	
		// Do a complete cycle
		for(int index = 0; index < J; index++) {
		
			double delta = ccdUpdateBeta(index);
			delta = applyBounds(delta, index);
			if (delta != 0) {
				sufficientStatisticsKnown = false;
				updateSufficientStatistics(delta, index);
			}
			
			if ((index+1) % 100 == 0) {
				//cout << "Finished variable " << (index+1) << endl;
			}
			
		}
		
		iteration++;
//		bool checkConvergence = (iteration % J == 0 || iteration == maxIterations);
		bool checkConvergence = true; // Check after each complete cycle

		if (checkConvergence) {

			double conv;
			if (convergenceType == LANGE) {
				double thisObjFunc = getObjectiveFunction();
				conv = computeConvergenceCriterion(thisObjFunc, lastObjFunc);
				lastObjFunc = thisObjFunc;
			} else { // ZHANG_OLES
				conv = computeZhangOlesConvergenceCriterion();
				saveXBeta();
			} // Necessary to call getObjFxn or computeZO before getLogLikelihood,
			  // since these copy over XBeta

			double thisLogLikelihood = getLogLikelihood();
			double thisLogPrior = getLogPrior();
			double thisLogPost = thisLogLikelihood + thisLogPrior;
			//cout << endl;
			//printVector(hBeta, J, cout);
			//cout << endl;
			//cout << "log post: " << thisLogPost
			//	 << " (" << thisLogLikelihood << " + " << thisLogPrior
			 //    << ") (iter:" << iteration << ") ";

			if (epsilon > 0 && conv < epsilon) {
			//	cout << "Reached convergence criterion" << endl;
				done = true;
			} else if (iteration == maxIterations) {
			//	cout << "Reached maximum iterations" << endl;
				done = true;
			} else {
			//	cout << endl;
			}
		}				
	}
}

/**
 * Computationally heavy functions
 */

void CyclicCoordinateDescent::computeGradientAndHession(int index, double *ogradient,
		double *ohessian) {
	bsccs::real gradient = 0;
	bsccs::real hessian = 0;

	int* nEvents = hNEvents;
	const int* end = hNEvents + N;

#ifdef MERGE_TRANSFORMATION
	bsccs::real* num = numerPid;
	bsccs::real* denom = denomPid;
	for (; nEvents != end; ++nEvents, ++num, ++denom) {
		const bsccs::real t = *num / *denom;
		const bsccs::real g = *nEvents * t;
		gradient += g;
		hessian += g * (static_cast<bsccs::real>(1.0) - t);
	}
#else
	bsccs::real* tmp = t1;
	for (; nEvents != end; ++nEvents, ++tmp) {
		const bsccs::real t = *tmp;
		const bsccs::real g = *nEvents * t;
		gradient += g;
		hessian += g * (static_cast<bsccs::real>(1.0) - t);
	}
#endif

	gradient -= hXjEta[index];
	*ogradient = static_cast<double>(gradient);
	*ohessian = static_cast<double>(hessian);
}

void CyclicCoordinateDescent::computeNumeratorForGradient(int index) {

	/*
	cout <<" GOT HERE " << endl;

	cout <<"Rows = " << hXI->getNumberOfRows() << endl;
	cout <<"Columns = " << hXI->getNumberOfColumns() << endl;
	cout <<hXI->getCompressedColumnVector(0)[2]<< endl;
	cout <<hXI->getCompressedColumnVector(1)[2]<< endl;
	cout <<hXI->getCompressedColumnVector(2)[2]<< endl;
	cout <<hXI->getCompressedColumnVector(3)[2]<< endl;
	cout <<hXI->getCompressedColumnVector(4)[2]<< endl;
*/
	zeroVector(numerPid, N);
	const int n = hXI->getNumberOfEntries(index);

#ifdef TEST_ROW_INDEX
	const int* rows = hXColumnRowIndicators[index];
#else
	const int* indicators = hXI->getCompressedColumnVector(index);
#endif

	for (int i = 0; i < n; i++) { // Loop through non-zero entries only
#ifdef TEST_ROW_INDEX
		const int k = rows[i];
		numerPid[k] += offsExpXBeta[k];
#else
		const int k = indicators[i];
		numerPid[hPid[k]] += offsExpXBeta[k];
#endif
	}
}

void CyclicCoordinateDescent::computeRatio(int index) {
	const int n = hXI->getNumberOfEntries(index);

#ifdef BETTER_LOOPS
	bsccs::real* t = t1;
	const bsccs::real* end = t + N;
	bsccs::real* num = numerPid;
	bsccs::real* denom = denomPid;
	for (; t != end; ++t, ++num, ++denom) {
		*t = *num / *denom;
	}
#else
	for (int i = 0; i < N; i++) {
		t1[i] = numerPid[i] / denomPid[i];
	}
#endif
}

void CyclicCoordinateDescent::computeRatiosForGradientAndHessian(int index) {
	computeNumeratorForGradient(index);
#ifndef MERGE_TRANSFORMATION
	computeRatio(index);
#endif
}

// tshaddox code...

double CyclicCoordinateDescent::getGradient(int drug) {
	double t1 = 1/sigma2Beta; // this is the hyperparameter that is used in the original code
	double t2 = 1/classHierarchyVariance;


	int parent = getParentMapCCD[drug];
	vector<int> siblings = getChildMapCCD[parent];
	double sumBetas = 0;
	int nSiblingsOfInterest = 0; //Different from siblings length if weights used
	//cout << "drug is " << drug << endl;
	//cout << "siblings are: " << endl;
	for (int i = 0; i < siblings.size(); i++) {
		//cout << siblings[i] << endl;
		sumBetas += hBeta[siblings[i]];
		//cout << "hBeta[" << siblings[i] << "] = " << hBeta[siblings[i]]<< endl;
	}


	double gradient = t1*hBeta[drug] - t1*t1*sumBetas / (siblings.size()*t1 + t2);
	//cout << "gradient = " << gradient << endl;
	return(gradient);
}

double CyclicCoordinateDescent::getHessian(int drug) {
	double t1 = 1/sigma2Beta;  // this is the hyperparameter used in the original code
	double t2 = 1/classHierarchyVariance;
	int parent = getParentMapCCD[drug];
	vector<int> siblings = getChildMapCCD[parent];
	double hessian = t1 - t1 / (siblings.size() + t2/t1);
	//cout << "hessian = " << hessian << endl;
	//cout << "1/sigma2Beta = " << 1.0/sigma2Beta << endl;
	return(hessian);
}

double CyclicCoordinateDescent::getGradient_lasso(int drug) {
	double groupLambda = sqrt(2.0 / classHierarchyVariance);
	double l2GroupSumSquared = 0;
	int parent = getParentMapCCD[drug];
	vector<int> siblings = getChildMapCCD[parent];
	for (int i = 0; i < siblings.size(); i++) {
		l2GroupSumSquared += hBeta[siblings[i]]*hBeta[siblings[i]];
	}
	double l2GroupNorm = sqrt(l2GroupSumSquared);
	double gradient;
	if (l2GroupNorm == 0) {
		gradient = groupLambda*1;//sign(hBeta[drug]);
	} else {
		gradient = groupLambda*(sqrt(hBeta[drug]*hBeta[drug])/l2GroupNorm); //groupLambda*(hBeta[drug]/l2GroupNorm);
	}
	return(gradient);
}

double CyclicCoordinateDescent::getHessian_lasso(int drug) {
	double groupLambda = sqrt(2.0 / classHierarchyVariance);
	double l2GroupSumSquared = 0;
	int parent = getParentMapCCD[drug];
	vector<int> siblings = getChildMapCCD[parent];
	for (int i = 0; i < siblings.size(); i++) {
		l2GroupSumSquared += hBeta[siblings[i]]*hBeta[siblings[i]];
	}
	double l2GroupNorm = sqrt(l2GroupSumSquared);
	double hessian;
	if (l2GroupNorm == 0) {
		hessian = 0;
	} else {
		hessian = groupLambda*(1/l2GroupNorm)*(1 - (hBeta[drug]*hBeta[drug]/l2GroupNorm));
	}
	return(hessian);
}

// tshaddox code end...

double CyclicCoordinateDescent::ccdUpdateBeta(int index) {

	double delta;

	if (!sufficientStatisticsKnown) {
		cerr << "Error in state synchronization." << endl;
		exit(0);		
	}
	
	computeRatiosForGradientAndHessian(index);
	
	double g_d1;
	double g_d2;
					
	computeGradientAndHession(index, &g_d1, &g_d2);
	if (priorType == NORMAL) {
		//cout << "NORMAL" << endl;
		if(useHierarchy) {
			delta = - (g_d1 + (getGradient(index))) /
					(g_d2 + getHessian(index));
		} else {
			delta = - (g_d1 + (hBeta[index] / sigma2Beta)) /
				  (g_d2 + (1.0 / sigma2Beta));
		}

		
	} else {
		//cout << "Laplacian" << endl;
		double neg_update;
		double pos_update;

		// tshaddox code (CHECK THIS)
		if(useHierarchy) {
			neg_update = - (g_d1 - lambda - getGradient_lasso(index)) /  (g_d2 + getHessian_lasso(index));
			pos_update = - (g_d1 + lambda + getGradient_lasso(index)) / (g_d2 + getHessian_lasso(index));
		} else {
			neg_update = - (g_d1 - lambda) / g_d2;
			pos_update = - (g_d1 + lambda) / g_d2;
		}
		int signBetaIndex = sign(hBeta[index]);
		
		if (signBetaIndex == 0) {
						
			if (neg_update < 0) {
				delta = neg_update;
			} else if (pos_update > 0) {
				delta = pos_update;
			} else {
				delta = 0;
			}
		} else { // non-zero entry
			
			if (signBetaIndex < 0) {
				delta = neg_update;
			} else {
				delta = pos_update;			
			}
			
			if ( sign(hBeta[index] + delta) != signBetaIndex ) {
				delta = - hBeta[index];
			}			
		}
	}
	//cout << "drug is " << index << " hbeta is " << hBeta[index] << " delta is " << delta << endl;
	return delta;
}

void CyclicCoordinateDescent::computeXBeta(void) {
	// Separate function for benchmarking
	cerr << "Computing X %*% beta" << endl;
	cerr << "Not yet implemented." << endl;
	exit(-1);
//	hXBeta = hX * hBeta;
}

void CyclicCoordinateDescent::updateXBeta(double delta, int index) {
	// Separate function for benchmarking
	hBeta[index] += delta;

	bsccs::real realDelta = static_cast<bsccs::real>(delta);


#ifdef TEST_ROW_INDEX
	const int* rows = hXColumnRowIndicators[index];
#endif

	const int* indicators = hXI->getCompressedColumnVector(index);
	const int n = hXI->getNumberOfEntries(index);
	for (int i = 0; i < n; i++) { // Loop through non-zero entries only
		const int k = indicators[i];
		hXBeta[k] += realDelta;
#ifdef TEST_SPARSE
		bsccs::real oldEntry = offsExpXBeta[k];
		bsccs::real newEntry = offsExpXBeta[k] = hOffs[k] * exp(hXBeta[k]);
#ifdef TEST_ROW_INDEX
		denomPid[rows[i]] += (newEntry - oldEntry);
#else
		denomPid[hPid[k]] += (newEntry - oldEntry);
#endif
#endif
	}	
}

void CyclicCoordinateDescent::updateSufficientStatistics(double delta, int index) {
	updateXBeta(delta, index);
	computeRemainingStatistics(false);
}

void CyclicCoordinateDescent::computeRemainingStatistics(bool allStats) {
	// Separate function for benchmarking
#ifdef TEST_SPARSE		
	if (allStats) {
#endif	
	zeroVector(denomPid, N);
	for (int i = 0; i < K; i++) {
		offsExpXBeta[i] = hOffs[i] * exp(hXBeta[i]);	
		denomPid[hPid[i]] += offsExpXBeta[i];
	}
#ifdef TEST_SPARSE			
	}
#endif	
	sufficientStatisticsKnown = true;
}

double CyclicCoordinateDescent::oneNorm(bsccs::real* vector, const int length) {
	double norm = 0;
	for (int i = 0; i < length; i++) {
		norm += abs(vector[i]);
	}
	return norm;
}

double CyclicCoordinateDescent::twoNormSquared(bsccs::real * vector, const int length) {
	double norm = 0;
	for (int i = 0; i < length; i++) {
		norm += vector[i] * vector[i];
	}
	return norm;
}

/**
 * Updating and convergence functions
 */

double CyclicCoordinateDescent::computeConvergenceCriterion(double newObjFxn, double oldObjFxn) {	
	// This is the stopping criterion that Ken Lange generally uses
	return abs(newObjFxn - oldObjFxn) / (abs(newObjFxn) + 1.0);
}

double CyclicCoordinateDescent::applyBounds(double inDelta, int index) {
	double delta = inDelta;
	if (delta < -hDelta[index]) {
		delta = -hDelta[index];
	} else if (delta > hDelta[index]) {
		delta = hDelta[index];
	}

	hDelta[index] = max(2.0 * abs(delta), 0.5 * hDelta[index]);
	return delta;
}

/**
 * Utility functions
 */

void CyclicCoordinateDescent::computeXjEta(void) {

	for (int drug = 0; drug < J; drug++) {
		hXjEta[drug] = 0;
		const int* indicators = hXI->getCompressedColumnVector(drug);
		const int n = hXI->getNumberOfEntries(drug);

		if (useCrossValidation) {
			for (int i = 0; i < n; i++) { // Loop through non-zero entries only
				const int k = indicators[i];
				hXjEta[drug] += hEta[k] * hWeights[k];
			}
		} else {
			for (int i = 0; i < n; i++) { // Loop through non-zero entries only
				const int k = indicators[i];
				hXjEta[drug] += hEta[k];
			}
		}
	}
}

template <class T>
void CyclicCoordinateDescent::printVector(T* vector, int length, ostream &os) {
	os << "(" << vector[0];
	for (int i = 1; i < length; i++) {
		os << ", " << vector[i];
	}
	os << ")";			
}

template <class T> 
T* CyclicCoordinateDescent::readVector(
		const char *fileName,
		int *length) {

	ifstream fin(fileName);
	T d;
	std::vector<T> v;
	
	while (fin >> d) {
		v.push_back(d);
	}
	
	T* ptr = (T*) malloc(sizeof(T) * v.size());
	for (int i = 0; i < v.size(); i++) {
		ptr[i] = v[i];
	}
		
	*length = v.size();
	return ptr;
}

void CyclicCoordinateDescent::testDimension(int givenValue, int trueValue, const char *parameterName) {
	if (givenValue != trueValue) {
		cerr << "Wrong dimension in " << parameterName << " vector." << endl;
		exit(-1);
	}
}	

inline int CyclicCoordinateDescent::sign(double x) {
	if (x == 0) {
		return 0;
	}
	if (x < 0) {
		return -1;
	}
	return 1;
}
}
