/*
 * ImputeVaraibles.h
 *
 *  Created on: Jul 28, 2012
 *      Author: Sushil Mittal
 */

#ifndef IMPUTEVARIABLES_H_
#define IMPUTEVARIABLES_H_

#include "ImputeInputReader.h"
#include "ccd.h"

class ImputeVariables {
public:
	ImputeVariables();
	~ImputeVariables();
	void initialize(CCDArguments arguments);
	void getComplement(vector<real>& weights);
	void impute(CCDArguments arguments);
	void randomizeImputationsLR(real* yPred, real* weights, int col);
	void randomizeImputationsLS(real* yPred, real* weights, int col);
private:
	CyclicCoordinateDescent*  ccd;
	AbstractModelSpecifics* model;
	InputReader* reader;
	int nImputations;
};


#endif