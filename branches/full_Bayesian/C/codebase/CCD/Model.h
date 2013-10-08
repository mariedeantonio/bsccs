/*
 * Model.h
 *
 *  Created on: Jul 24, 2013
 *      Author: tshaddox
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include "Parameter.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "CyclicCoordinateDescent.h"


namespace bsccs {
class Model {
public:
	Model();
	virtual ~Model();

	void generateCholesky();
	void initialize(CyclicCoordinateDescent& ccd);
	void initializeHessian();
	void clearHessian();

	void Beta_HatRestorableSet(bool restorable);
	void BetaRestorableSet(bool restorable);
	void SigmaSquaredRestorableSet(bool restorable);

	void restore();

	void logState();

	void resetWithNewSigma(CyclicCoordinateDescent& ccd);

	bool getUseHastingsRatio();
	void setUseHastingsRatio(bool useHastingsRatio);

	double getLoglikelihood();
	void setLoglikelihood(double useHastingsRatio);



	Parameter& getBeta();
	Parameter& getBeta_Hat();
	Parameter& getSigmaSquared();
	Eigen::LLT<Eigen::MatrixXf> getCholeskyLLT();
	Eigen::MatrixXf& getHessian();

private:
	Eigen::LLT<Eigen::MatrixXf> CholDecom;

	Eigen::MatrixXf HessianMatrix;

	vector<vector<bsccs::real> > hessian;

	Parameter Beta_Hat;
	Parameter Beta;
	Parameter SigmaSquared;

	double loglikelihood;

	bool Beta_HatRestorable;
	bool BetaRestorable;
	bool SigmaSquaredRestorable;
	bool useHastingsRatio;

	int J;


};
}





#endif /* MODEL_H_ */
