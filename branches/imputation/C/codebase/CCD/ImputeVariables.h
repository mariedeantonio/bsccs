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
private:
	CyclicCoordinateDescent*  ccd;
	AbstractModelSpecifics* model;
	ImputeInputReader* reader;
	int nImputations;
};


#endif