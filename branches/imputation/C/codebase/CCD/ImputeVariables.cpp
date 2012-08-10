/*
 * ImputeVaraibles.cpp
 *
 *  Created on: Jul 28, 2012
 *      Author: Sushil Mittal
 */

#include "ImputeVariables.h"
#include <time.h>

double randn_notrig(double mu=0.0, double sigma=1.0) {
	static bool deviateAvailable=false;        //        flag
	static float storedDeviate;                        //        deviate from previous calculation
	double polar, rsquared, var1, var2;

	//        If no deviate has been stored, the polar Box-Muller transformation is
	//        performed, producing two independent normally-distributed random
	//        deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable) {

		//        choose pairs of uniformly distributed deviates, discarding those
		//        that don't fall within the unit circle
		do {
			var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			rsquared=var1*var1+var2*var2;
		} while ( rsquared>=1.0 || rsquared == 0.0);

		//        calculate polar tranformation for each deviate
		polar=sqrt(-2.0*log(rsquared)/rsquared);

		//        store first deviate and set flag
		storedDeviate=var1*polar;
		deviateAvailable=true;

		//        return second deviate
		return var2*polar*sigma + mu;
	}

	//        If a deviate is available from a previous call to this function, it is
	//        returned, and the flag is set to false.
	else {
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
}


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
	
	srand(time(NULL));
	
	vector<int> nMissingPerColumn = reader->getnMissingPerColumn();
	vector<string> columnTypes = reader->getColumnTypesToImpute();
	
	for(int i = 0; i < nImputations; i++){
		for(int j = 0; j < (int)nMissingPerColumn.size(); j++){
			if(nMissingPerColumn[j] > 0){
				vector<real> weights;
				reader->setupDataForImputation(j,weights);
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
				
				int nRows = reader->getNumberOfRows();
				vector<real> yPred(nRows);
				vector<real> allOnes(nRows,1.0);

				ccd->getPredictiveEstimates(&yPred[0], &allOnes[0]);
				
				if(columnTypes[j]=="ls")
					randomizeImputationsLS(&yPred[0], &weights[0], j);
				else
					randomizeImputationsLR(&yPred[0], &weights[0], j);

				if (ccd)
					delete ccd;
				if (model)
					delete model;
			}
		}
	}
}

void ImputeVariables::randomizeImputationsLR(real* yPred, real* weights, int col){
	real* y = reader->getDataVector(col);
	for(int i = 0; i < reader->getNumberOfRows(); i++){
		if(weights[i]){
			if(yPred[i] < rand()/RAND_MAX)
				y[i] = 1.0;
			else
				y[i] = 0.0;
		}
	}
}

void ImputeVariables::randomizeImputationsLS(real* yPred, real* weights, int col){
	
	int n = 0;
	vector<real> xMean(col);
	vector<real> xVar(col);
	
	int nRows = reader->getNumberOfRows();
	vector<real> dist(nRows,0.0);

	reader->getSampleMeanVariance(col,&xMean[0],&xVar[0]);
	real* y = reader->getDataVector(col);

	real sigma = 0.0;
	for(int i = 0; i < nRows; i++){
		if(!weights[i]){
			sigma += (y[i] - yPred[i]) * (y[i] - yPred[i]);
			n++;
		}
		else{
			vector<real> x(col);
			reader->getDataRow(i,&x[0]);
			for(int j = 0; j < col; j++){
				if(xVar[j] > 1e-10)
					dist[i] += (x[j] - xMean[j]) * (x[j] - xMean[j]) / xVar[j] ;
			}
		}
	}
	sigma = sigma/(n-2);

	for(int i = 0; i < nRows; i++){
		if(weights[i]){
			double predInterval = 1.96 * sigma * sqrt(1 + 1/n + dist[i]/(n-1));
			y[i] = randn_notrig(yPred[i],predInterval);
		}
	}
}