/*
 * ImputeVaraibles.cpp
 *
 *  Created on: Jul 28, 2012
 *      Author: Sushil Mittal
 */

#include "ImputeVariables.h"


ImputeVariables::ImputeVariables(){
	ccd = NULL;
	model = NULL;
}

ImputeVariables::~ImputeVariables(){
	if (reader)
		delete reader;
}

void ImputeVariables::initialize(CCDArguments arguments){
	nImputations = arguments.numberOfImputations;
	reader = new ImputeInputReader();
	reader->readFile(arguments.inFileName.c_str());
}

void ImputeVariables::getComplement(vector<real>& weights){
	for(std::vector<real>::iterator it = weights.begin(); it != weights.end(); it++) {
		*it = 1 - *it;
	}
}

void ImputeVariables::impute(CCDArguments arguments){
	
	vector<int> nMissingPerColumn = reader->getnMissingPerColumn();
	vector<string> columnTypes = reader->getColumnTypesToImpute();
	for(int i = 0; i < nImputations; i++){
		for(int j = 0; j < (int)nMissingPerColumn.size(); j++){
			if(nMissingPerColumn[j] > 0){
				vector<real> weights;
				reader->setupData(j,weights);
				if(columnTypes[j]=="ls")	//Least Squares
					model = new ModelSpecifics<LeastSquares<real>,real>();
				else						//Logistic Regression
					model = new ModelSpecifics<LogisticRegression<real>,real>();

				ccd = new CyclicCoordinateDescent(reader, *model);
				ccd->setWeights(&weights[0]);

				if(arguments.useNormalPrior)
					ccd->setPriorType(NORMAL);
				if (arguments.hyperPriorSet)
					ccd->setHyperprior(arguments.hyperprior);

				ccd->update(arguments.maxIterations, arguments.convergenceType, arguments.tolerance);
				getComplement(weights);
				vector<real> y = ccd->getPredictiveEstimates(&weights[0]);

				//TODO: fill data[col]->at(absententries) with y
				if (ccd)
					delete ccd;
				if (model)
					delete model;
			}
		}
	}
}