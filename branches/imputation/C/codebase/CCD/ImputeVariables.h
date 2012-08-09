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
	void randomizeImputationsLR(real* y, real* weights, int n);
	void randomizeImputationsLS(real* y, real* weights, int n);
private:
	CyclicCoordinateDescent*  ccd;
	AbstractModelSpecifics* model;
	ImputeInputReader* reader;
	int nImputations;
};


#endif