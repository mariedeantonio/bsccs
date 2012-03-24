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
#include <set>

#include "CyclicCoordinateDescent.h"
#include "InputReader.h"
#include "Iterators.h"

#define PI	3.14159265358979323851280895940618620443274267017841339111328125

using namespace std;

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

	hXI = reader;
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
//	free(t1);
	
#ifdef NO_FUSE
	free(wPid);
#endif
	
	if (hWeights) {
		free(hWeights);
	}

#ifdef SPARSE_PRODUCT
	for (std::vector<std::vector<int>* >::iterator it = sparseIndices.begin(); it != sparseIndices.end(); ++it) {
		if (*it) {
			delete *it;
		}
	}
#endif
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
	hDelta = (real*) malloc(J * sizeof(real));
//	for (int j = 0; j < J; j++) {
//		hDelta[j] = 2.0;
//	}

	hBeta = (real*) calloc(J, sizeof(real)); // Fixed starting state
	hXBeta = (real*) calloc(K, sizeof(real));
	hXBetaSave = (real*) calloc(K, sizeof(real));

	// Set prior
	priorType = LAPLACE;
	sigma2Beta = 1000;
	lambda = sqrt(2.0/20.0);
	
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
	offsExpXBeta = (real*) malloc(sizeof(real) * K);
	xOffsExpXBeta = (real*) malloc(sizeof(real) * K);

	// Put numer and denom in single memory block, with first entries on 16-word boundary
	int alignedLength = getAlignedLength(N);

	numerPid = (real*) malloc(sizeof(real) * 3 * alignedLength);
//	denomPid = (real*) malloc(sizeof(real) * N);
	denomPid = numerPid + alignedLength; // Nested in denomPid allocation
	numerPid2 = numerPid + 2 * alignedLength;
//	t1 = (real*) malloc(sizeof(real) * N);
	hNEvents = (int*) malloc(sizeof(int) * N);
	hXjEta = (real*) malloc(sizeof(real) * J);
	hWeights = NULL;
	
#ifdef NO_FUSE
	wPid = (real*) malloc(sizeof(real) * alignedLength);
#endif

	for (int j = 0; j < J; ++j) {
		if (hXI->getFormatType(j) == DENSE) {
			sparseIndices.push_back(NULL);
		} else {
			std::set<int> unique;
			const int n = hXI->getNumberOfEntries(j);
			const int* indicators = hXI->getCompressedColumnVector(j);
			for (int j = 0; j < n; j++) { // Loop through non-zero entries only
				const int k = indicators[j];
				const int i = hPid[k];
				unique.insert(i);
			}
			std::vector<int>* indices = new std::vector<int>(unique.begin(),
					unique.end());
			sparseIndices.push_back(indices);
		}
	}

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

	outLog << "Drug_concept_id" << sep << "Condition_concept_id" << sep << "score" << endl;

	for (int i = 0; i < J; i++) {
		outLog << drugMap[i] << sep <<
//		i << sep <<
		conditionId << sep << hBeta[i] << endl;
	}
	outLog.close();
}

double CyclicCoordinateDescent::getPredictiveLogLikelihood(real* weights) {

	if (!sufficientStatisticsKnown) {
		computeRemainingStatistics(true, 0); // TODO Remove index????
	}

	getDenominators();

	double logLikelihood = 0;

	for (int i = 0; i < K; i++) {
		logLikelihood += hEta[i] * weights[i] * (hXBeta[i] - log(denomPid[hPid[i]]));
	}

	return logLikelihood;
}

int CyclicCoordinateDescent::getBetaSize(void) {
	return J;
}

real CyclicCoordinateDescent::getBeta(int i) {
	if (!sufficientStatisticsKnown) {
		computeRemainingStatistics(true, i);
	}
	return hBeta[i];
}

double CyclicCoordinateDescent::getLogLikelihood(void) {

	if (!sufficientStatisticsKnown) {
		computeRemainingStatistics(true, 0); // TODO Check index?
	}

	getDenominators();

	real logLikelihood = 0;

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
		logLikelihood -= hNEvents[i] * log(denomPid[i]); // TODO What about weights?
	}

	return static_cast<double>(logLikelihood);
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

void CyclicCoordinateDescent::setWeights(real* iWeights) {

	if (iWeights == NULL) {
		std::cerr << "Turning off weights!" << std::endl;
		// Turn off weights
		useCrossValidation = false;
		validWeights = false;
		sufficientStatisticsKnown = false;
		return;
	}

	if (hWeights == NULL) {
		hWeights = (real*) malloc(sizeof(real) * K);
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
	real criterion = 0;
	if (useCrossValidation) {
		for (int i = 0; i < K; i++) {
			criterion += hXBeta[i] * hEta[i] * hWeights[i];
		}
	} else {
		for (int i = 0; i < K; i++) {
			criterion += hXBeta[i] * hEta[i];
		}
	}
	return static_cast<double>(criterion);
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
	memcpy(hXBetaSave, hXBeta, K * sizeof(real));
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
		computeRemainingStatistics(true, 0); // TODO Check index?
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
			if (delta != 0.0) {
				sufficientStatisticsKnown = false;
				updateSufficientStatistics(delta, index);
			}
			
			if ((index+1) % 100 == 0) {
				cout << "Finished variable " << (index+1) << endl;
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
			cout << endl;
			printVector(hBeta, J, cout);
			cout << endl;
			cout << "log post: " << thisLogPost
				 << " (" << thisLogLikelihood << " + " << thisLogPrior
			     << ") (iter:" << iteration << ") ";

			if (epsilon > 0 && conv < epsilon) {
				cout << "Reached convergence criterion" << endl;
				done = true;
			} else if (iteration == maxIterations) {
				cout << "Reached maximum iterations" << endl;
				done = true;
			} else {
				cout << endl;
			}
		}				
	}
}

/**
 * Computationally heavy functions
 */

void CyclicCoordinateDescent::computeGradientAndHessian(int index, double *ogradient,
		double *ohessian) {
	// Run-time dispatch
	switch (hXI->getFormatType(index)) {
	case INDICATOR :
		computeGradientAndHessianImpl<IndicatorIterator>(index, ogradient, ohessian);
		break;
	case SPARSE :
		computeGradientAndHessianImpl<SparseIterator>(index, ogradient, ohessian);
		break;
	case DENSE :
//		computeGradientAndHessianImplHand(index, ogradient, ohessian); // TODO This appears slightly faster
		computeGradientAndHessianImpl<DenseIterator>(index, ogradient, ohessian);
		break;
	}
}

//void CyclicCoordinateDescent::computeGradientAndHessianImplHand(int index, double *ogradient,
//		double *ohessian) {
//	real gradient = 0;
//	real hessian = 0;
//	for (int k = 0; k < N; ++k) {
//		const real t = numerPid[k] / denomPid[k];
//		const real g = hNEvents[k] * t;
//		gradient += g;
//		hessian += g * (static_cast<real>(1.0) - t); // TODO Update for !indicators
//	}
//
//	gradient -= hXjEta[index];
//	*ogradient = static_cast<double>(gradient);
//	*ohessian = static_cast<double>(hessian);
//}

template <>
inline void CyclicCoordinateDescent::incrementGradientAndHessian<IndicatorIterator>(
		real* gradient, real* hessian,
		real numer, real numer2, real denom, int nEvents) {
	const real t = numer / denom;
	const real g = nEvents * t;
	*gradient += g;
	*hessian += g * (static_cast<real>(1.0) - t);
}

template <class IteratorType>
inline void CyclicCoordinateDescent::incrementGradientAndHessian(
		real* gradient, real* hessian,
		real numer, real numer2, real denom, int nEvents) {
	const real t = numer / denom;
	const real g = nEvents * t;
	*gradient += g;
	if (IteratorType::isIndicator) {
		*hessian += g * (static_cast<real>(1.0) - t);
	} else {
//		const real h1 = nEvents * (numer2 * denom - numer * numer) / (denom * denom);
		const real h1 = nEvents * (numer2 / denom - t * t); // Bounded by x_j^2
		const real h2 =  g * (static_cast<real>(1.0) - t);
//		cerr << "Not indy! " << h1 << " " << h2 << endl;
//		exit(-1);
		*hessian += h1;
	}
}

//template <>
//inline void CyclicCoordinateDescent::incrementGradientAndHessian<SparseIterator>(
//		real* gradient, real* hessian,
//		real numer, real numer2, real denom, int nEvents) {
//	const real t = numer / denom;
//	const real g = nEvents * t;
//	*gradient += g;
////	const real h1 = nEvents * (numer2 * denom - numer * numer) / (denom * denom);
//	const real h1 = nEvents * (numer2 / denom - t * t);
//	const real h2 =  g * (static_cast<real>(1.0) - t);
//	*hessian += g * (static_cast<real>(1.0) - t);
//	cerr << "Sparse it! " << h1 << " " << h2 << endl;
//	exit(-1);
//}
//
//template <>
//inline void CyclicCoordinateDescent::incrementGradientAndHessian<DenseIterator>(
//		real* gradient, real* hessian,
//		real numer, real numer2, real denom, int nEvents) {
//	const real t = numer / denom;
//	const real g = nEvents * t;
//	*gradient += g;
//	const real h1 = nEvents * (numer2 * denom - numer * numer) / (denom * denom);
//	const real h2 =  g * (static_cast<real>(1.0) - t);
//	*hessian += g * (static_cast<real>(1.0) - t);
//	cerr << "Dense it! " << h1 << " " << h2 << endl;
////	exit(-1);
//}

template <class IteratorType>
void CyclicCoordinateDescent::computeGradientAndHessianImpl(int index, double *ogradient,
		double *ohessian) {
	real gradient = 0;
	real hessian = 0;
	
	IteratorType it(*sparseIndices[index], N); // TODO How to create with different constructor signatures?
	for (; it; ++it) {
		const int k = it.index();
		incrementGradientAndHessian<IteratorType>(
				&gradient, &hessian,
				numerPid[k], numerPid2[k], denomPid[k], hNEvents[k]
		);
	}
	gradient -= hXjEta[index];
	
	*ogradient = static_cast<double>(gradient);
	*ohessian = static_cast<double>(hessian);
}

//#define MD

void CyclicCoordinateDescent::computeNumeratorForGradient(int index) {
	// Run-time delegation
	switch (hXI->getFormatType(index)) {
		case INDICATOR : {
			IndicatorIterator it(*sparseIndices[index]);
			for (; it; ++it) { // Only affected entries
				numerPid[it.index()] = static_cast<real>(0.0);
			}
			incrementNumeratorForGradientImpl<IndicatorIterator>(index); }
			break;
		case DENSE :
			zeroVector(numerPid, N);
			zeroVector(numerPid2, N);
			incrementNumeratorForGradientImpl<DenseIterator>(index);
			break;
		case SPARSE : {
			IndicatorIterator it(*sparseIndices[index]);
			for (; it; ++it) { // Only affected entries
				numerPid[it.index()] = static_cast<real>(0.0);
				numerPid2[it.index()] = static_cast<real>(0.0); // TODO Does this invalid the cache line too much?
			}
			incrementNumeratorForGradientImpl<SparseIterator>(index); }
			break;
		default :
			// throw error
			exit(-1);
	}
}

template <class IteratorType>
void CyclicCoordinateDescent::incrementNumeratorForGradientImpl(int index) {
	IteratorType it(*hXI, index);
	for (; it; ++it) {
		const int k = it.index();
		numerPid[hPid[k]] += offsExpXBeta[k] * it.value();
		if (!IteratorType::isIndicator) { // TODO Ensure compile-time evaluation
//			cerr << "Here?" << endl;
//			exit(-1);
			numerPid2[hPid[k]] += offsExpXBeta[k] * it.value() * it.value();
//			cerr << "numer2 = " << numerPid2[hPid[k]] << endl;
		}
	}
}

//void CyclicCoordinateDescent::computeRatio(int index) {
//	const int n = hXI->getNumberOfEntries(index);
//
//#ifdef BETTER_LOOPS
//	real* t = t1;
//	const real* end = t + N;
//	real* num = numerPid;
//	real* denom = denomPid;
//	for (; t != end; ++t, ++num, ++denom) {
//		*t = *num / *denom;
//	}
//#else
//	for (int i = 0; i < N; i++) {
//		t1[i] = numerPid[i] / denomPid[i];
//	}
//#endif
//}

void CyclicCoordinateDescent::computeRatiosForGradientAndHessian(int index) {
//	computeNumeratorForGradient(index);
//#ifndef MERGE_TRANSFORMATION
//	computeRatio(index);
//#endif
	cerr << "Error!" << endl;
	exit(-1);
}

double CyclicCoordinateDescent::ccdUpdateBeta(int index) {

	double delta;

	if (!sufficientStatisticsKnown) {
		cerr << "Error in state synchronization." << endl;
		exit(0);		
	}
	
	computeNumeratorForGradient(index);
	
	double g_d1;
	double g_d2;
					
	computeGradientAndHessian(index, &g_d1, &g_d2);

//	cerr << "grad = " << g_d1 << endl;
//	cerr << "hess = " << g_d2 << endl;
//	exit(-1);

	// TODO Remove check
	if (g_d2 <= 0.0) {
		cerr << "Not positive-definite!" << endl;
		exit(-1);
	}
		
	// Move into separate delegate-function (below)
	if (priorType == NORMAL) {
		
		delta = - (g_d1 + (hBeta[index] / sigma2Beta)) /
				  (g_d2 + (1.0 / sigma2Beta));
		
	} else {
					
		double neg_update = - (g_d1 - lambda) / g_d2;
		double pos_update = - (g_d1 + lambda) / g_d2;
		
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
	
	return delta;
}

void CyclicCoordinateDescent::computeXBeta(void) {
	// Separate function for benchmarking
	cerr << "Computing X %*% beta" << endl;
	cerr << "Not yet implemented." << endl;
	exit(-1);
//	hXBeta = hX * hBeta;
}

template <class IteratorType>
void CyclicCoordinateDescent::updateXBetaImpl(real realDelta, int index) {
	IteratorType it(*hXI, index);
	for (; it; ++it) {
		const int k = it.index();
		hXBeta[k] += realDelta * it.value();
		// Update denominators as well
		real oldEntry = offsExpXBeta[k];
		real newEntry = offsExpXBeta[k] = hOffs[k] * exp(hXBeta[k]);
		denomPid[hPid[k]] += (newEntry - oldEntry);
	}
}

void CyclicCoordinateDescent::updateXBetaImplHand(real realDelta, int index) {

	real* data = hXI->getDataVector(index);
	real* xBeta = hXBeta;
	real* offsEXB = offsExpXBeta;
	int* offs = hOffs;
	int* pid = hPid;

	for (int k = 0; k < K; k++) {
		*xBeta += realDelta * *data;
		real oldEntry = *offsEXB;
		real newEntry = *offsEXB = *offs * exp(*xBeta);
		denomPid[*pid] += (newEntry - oldEntry);
		data++; xBeta++; offsEXB++; offs++; pid++;
	}
}

void CyclicCoordinateDescent::updateXBeta(double delta, int index) {
	real realDelta = static_cast<real>(delta);
	hBeta[index] += realDelta;

	// Run-time dispatch to implementation depending on covariate FormatType
	switch(hXI->getFormatType(index)) {
		case INDICATOR :
			updateXBetaImpl<IndicatorIterator>(realDelta, index);
			break;
		case DENSE :
			updateXBetaImpl<DenseIterator>(realDelta, index);
//			updateXBetaImplHand(realDelta, index); // Oddly slower!
			break;
		case SPARSE :
			updateXBetaImpl<SparseIterator>(realDelta, index);
			break;
		default :
			// throw error
			exit(-1);
	}
}

void CyclicCoordinateDescent::updateSufficientStatistics(double delta, int index) {
	updateXBeta(delta, index);
	sufficientStatisticsKnown = true;
//	computeRemainingStatistics(false);
}

void CyclicCoordinateDescent::computeRemainingStatistics(bool allStats, int index) { // TODO Rename
	// Separate function for benchmarking

	if (allStats) {
		zeroVector(denomPid, N);
		for (int i = 0; i < K; i++) {
			offsExpXBeta[i] = hOffs[i] * exp(hXBeta[i]);
			denomPid[hPid[i]] += offsExpXBeta[i];
		}
//		cerr << "den[0] = " << denomPid[0] << endl;
//		exit(-1);
	}
	sufficientStatisticsKnown = true;
}

double CyclicCoordinateDescent::oneNorm(real* vector, const int length) {
	double norm = 0;
	for (int i = 0; i < length; i++) {
		norm += abs(vector[i]);
	}
	return norm;
}

double CyclicCoordinateDescent::twoNormSquared(real * vector, const int length) {
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

//	cerr << "YXj";
	for (int drug = 0; drug < J; drug++) {
		hXjEta[drug] = 0;
		GenericIterator it(*hXI, drug);

		if (useCrossValidation) {
			for (; it; ++it) {
				const int k = it.index();
				hXjEta[drug] += it.value() * hEta[k] * hWeights[k];
			}
		} else {
			for (; it; ++it) {
				const int k = it.index();
				hXjEta[drug] += it.value() * hEta[k];
			}
		}
//		cerr << " " << hXjEta[drug];
	}
//	cerr << endl;
//	exit(1);
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
