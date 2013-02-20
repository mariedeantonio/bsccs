/*
 * MCMCDriver.h
 *
 *  Created on: Aug 6, 2012
 *      Author: trevorshaddox
 */

#ifndef MCMCDRIVER_H_
#define MCMCDRIVER_H_


#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include "CyclicCoordinateDescent.h"
#include "CrossValidationSelector.h"
#include "AbstractSelector.h"
#include "ccd.h"

namespace bsccs {
class MCMCDriver {
public:
	MCMCDriver(InputReader * inReader);

	virtual ~MCMCDriver();

	virtual void drive(
			CyclicCoordinateDescent& ccd, double betaAmount);

	void generateCholesky();

	void initializeHessian();

private:
	InputReader* reader;

	double acceptanceTuningParameter;

	Eigen::LLT<Eigen::MatrixXf> CholDecom;

	double acceptanceRatioTarget;

	vector<vector<double> > MCMCResults_BetaVectors;

	void adaptiveKernel(int numberIterations, int numberAcceptances);

	vector<double> MCMCResults_SigmaSquared;

	vector<double> BetaValues;

	vector<vector<bsccs::real> > hessian;

	vector<vector<bsccs::real> > cholesky;

	int J;

	int maxIterations;

	int nBetaSamples;

	int nSigmaSquaredSamples;
};
}


#endif /* MCMCDRIVER_H_ */
