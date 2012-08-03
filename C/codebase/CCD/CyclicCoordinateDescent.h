/*
 * CyclicCoordinateDescent.h
 *
 *  Created on: May-June, 2010
 *      Author: msuchard
 */

#ifndef CYCLICCOORDINATEDESCENT_H_
#define CYCLICCOORDINATEDESCENT_H_

#include <Eigen/Dense>
#include "CompressedDataMatrix.h"
#include "InputReader.h"




//using namespace std;
//using namespace Eigen;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::ofstream;

//#define DEBUG

#define TEST_SPARSE // New sparse updates are great
//#define TEST_ROW_INDEX
#define BETTER_LOOPS
#define MERGE_TRANSFORMATION
#define NEW_NUMERATOR
#define SPARSE_PRODUCT

#define USE_ITER


//#define NO_FUSE

#ifdef DOUBLE_PRECISION
	typedef double realTRS;
#else
	typedef float realTRS;
#endif 


enum PriorType {
	LAPLACE = 0,
	NORMAL  = 1
};

enum ConvergenceType {
	LANGE = 0,
	ZHANG_OLES = 1
};

class CyclicCoordinateDescent {
	
public:
	
	CyclicCoordinateDescent(void);
	
	CyclicCoordinateDescent(			
			const char* fileNameX,
			const char* fileNameEta,
			const char* fileNameOffs,
			const char* fileNameNEvents,
			const char* fileNamePid			
		);
	
	CyclicCoordinateDescent(
			InputReader* reader
		);

	CyclicCoordinateDescent(
			int inN,
			CompressedDataMatrix* inX,
			int* inEta, 
			int* inOffs, 
			int* inNEvents,
			int* inPid
		);
	
	void logResults(const char* fileName);

	virtual ~CyclicCoordinateDescent();
	
	double getLogLikelihood(void);

	double getPredictiveLogLikelihood(realTRS* weights);

	double getLogPrior(void);
	
	virtual double getObjectiveFunction(void);

	realTRS getBeta(int i);

	int getBetaSize(void);
		
	void update(int maxIterations, int convergenceType, double epsilon);

	virtual void resetBeta(void);

	// Setters
	void setHyperprior(double value);

	void setPriorType(int priorType);

	void setWeights(realTRS* weights);

	void setLogisticRegression(bool idoLR);

	void getHessianForCholesky_GSL();

	void getHessianForCholesky_Eigen(Eigen::MatrixXf* ReturnHessian);

	void getCholeskyFromHessian(Eigen::MatrixXf* HessianMatrix, Eigen::MatrixXf* CholeskyDecomp);

//	template <typename T>
	void setBeta(const std::vector<double>& beta);

//	void double getHessianComponent(int i, int j);

	// Getters
	string getPriorInfo();

	string getConditionId() const {
		return conditionId;
	}

	int getUpdateCount() const {
		return updateCount;
	}

	int getLikelihoodCount() const {
		return likelihoodCount;
	}
		
protected:
	
//private:
	
	void init(void);
	
	void resetBounds(void);

	void computeXBeta(void);

	void saveXBeta(void);

	void computeXjEta(void);

	void computeSufficientStatistics(void);

	void updateSufficientStatistics(double delta, int index);

	template <class IteratorType>
	void incrementNumeratorForGradientImpl(int index);

	void incrementNumeratorForGradientImplHand(int index);

	void computeNumeratorForGradient(int index);

	virtual void computeNEvents(void);

	virtual void updateXBeta(double delta, int index);

	template <class IteratorType>
	void updateXBetaImpl(realTRS delta, int index);

	void updateXBetaImplHand(realTRS realDelta, int index);

	virtual void computeRemainingStatistics(bool skip, int index);
	
	virtual void computeRatiosForGradientAndHessian(int index);

//	virtual void computeRatio(int index);

	virtual void computeGradientAndHessian(
			int index,
			double *gradient,
			double *hessian);

	template <class IteratorType>
	void computeGradientAndHessianImpl(
			int index,
			double *gradient,
			double *hessian);

	void computeGradientAndHessianImplHand(
			int index,
						double *gradient,
						double *hessian);

	template <class IteratorType>
	inline realTRS computeHessian(
			realTRS numer, realTRS numer2, realTRS denom,
			realTRS g, realTRS t);

	template <class IteratorType>
	inline void incrementGradientAndHessian(
			realTRS* gradient, realTRS* hessian,
			realTRS numer, realTRS numer2, realTRS denom, int nEvents);

	template <class IteratorType>
	void axpy(realTRS* y, const realTRS alpha, const int index);

	virtual void getDenominators(void);

	double computeLogLikelihood(void);

	double ccdUpdateBeta(int index);
	
	double applyBounds(
			double inDelta,
			int index);
	
	double computeConvergenceCriterion(double newObjFxn, double oldObjFxn);
	
	virtual double computeZhangOlesConvergenceCriterion(void);

	template <class T>
	void fillVector(T* vector, const int length, const T& value) {
		for (int i = 0; i < length; i++) {
			vector[i] = value;
		}
	}

	template <class T>
	void zeroVector(T* vector, const int length) {
		for (int i = 0; i < length; i++) {
			vector[i] = 0;
		}
	}

	int getAlignedLength(int N);
		
	void testDimension(int givenValue, int trueValue, const char *parameterName);
	
	template <class T>
	void printVector(T* vector, const int length, ostream &os);
	
	double oneNorm(realTRS* vector, const int length);
	
	double twoNormSquared(realTRS * vector, const int length);
	
	int sign(double x); 
	
	template <class T> 
	T* readVector(const char *fileName, int *length); 
			
	// Local variables
	
	//InputReader* hReader;
	
	ofstream outLog;
	bool hasLog;

	CompressedDataMatrix* hXI; // K-by-J-indicator matrix

	int* hOffs;  // K-vector
	int* hEta; // K-vector
	int* hNEvents; // K-vector
	int* hPid; // N-vector
	int** hXColumnRowIndicators; // J-vector
 	
	realTRS* hBeta;
	realTRS* hXBeta;
	realTRS* hXBetaSave;
	realTRS* hDelta;

	int N; // Number of patients
	int K; // Number of exposure levels
	int J; // Number of drugs
	
	string conditionId;

	int priorType;
	double sigma2Beta;
	double lambda;

	realTRS denomNullValue;

	bool sufficientStatisticsKnown;
	bool xBetaKnown;

	bool validWeights;
	bool useCrossValidation;
	bool doLogisticRegression;
	realTRS* hWeights;

	// temporary variables
	realTRS* expXBeta;
	realTRS* offsExpXBeta;
	realTRS* denomPid;
	realTRS* numerPid;
	realTRS* numerPid2;
	realTRS* xOffsExpXBeta;
	realTRS* hXjEta;

	int updateCount;
	int likelihoodCount;

#ifdef SPARSE_PRODUCT
	std::vector<std::vector<int>* > sparseIndices;
#endif
	
#ifdef NO_FUSE
	realTRS* wPid;
#endif
};

double convertVarianceToHyperparameter(double variance);

#endif /* CYCLICCOORDINATEDESCENT_H_ */
