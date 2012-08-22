/*
 * CrossValidationDriver.cpp
 *
 *  Created on: Sep 10, 2010
 *      Author: msuchard
 */

#include <iostream>
#include <iomanip>
#include <numeric>
#include <math.h>
#include <cstdlib>

#include "CrossValidationDriver.h"

namespace bsccs {

CrossValidationDriver::CrossValidationDriver(
			int iGridSize,
			double iLowerLimit,
			double iUpperLimit) : gridSize(iGridSize),
			lowerLimit(iLowerLimit), upperLimit(iUpperLimit) {

	// Do anything???

}

CrossValidationDriver::~CrossValidationDriver() {
	// Do nothing
}

double CrossValidationDriver::computeGridPoint(int step) {
	if (gridSize == 1) {
		return upperLimit;
	}
	// Linear grid
//	double stepSize = (upperLimit - lowerLimit) / (gridSize - 1);
//	return lowerLimit + step * stepSize;
	// Log uniform grid
	double stepSize = (log(upperLimit) - log(lowerLimit)) / (gridSize - 1);
	return exp(log(lowerLimit) + step * stepSize);
}
double CrossValidationDriver::computePointEstimate(const std::vector<double>& value) {
	// Mean of log values
	return accumulate(value.begin(), value.end(), 0.0);
}


void CrossValidationDriver::logResults(const CCDArguments& arguments) {

	ofstream outLog(arguments.cvFileName.c_str());
	if (!outLog) {
		cerr << "Unable to open log file: " << arguments.cvFileName << endl;
		exit(-1);
	}

	string sep(","); // TODO Make option

	double maxPoint;
	double maxValue;
	double minValue;
	findMax(&maxPoint, &maxValue, &minValue);

	for (int i = 0; i < gridPoint.size(); i++) {
		outLog << std::setw(5) << std::setprecision(4) << std::fixed << gridPoint[i] << sep;
		if (!arguments.useNormalPrior) {
			outLog << convertVarianceToHyperparameter(gridPoint[i]) << sep;
		}
		outLog << std::scientific << gridValue[i] << sep;
		outLog << (maxValue - gridValue[i]) << std::endl;
	}

	outLog.close();
}

void CrossValidationDriver::resetForOptimal(
		CyclicCoordinateDescent& ccd,
		CrossValidationSelector& selector,
		const CCDArguments& arguments) {

	ccd.setWeights(NULL);

	double maxPoint;
	double maxValue;
	double minValue;
	findMax(&maxPoint, &maxValue, &minValue);
	ccd.setHyperprior(maxPoint);
	ccd.resetBeta(); // Cold-start
}

/* Author: tshaddox
 * Helper function to translate effects in the greedy cross validation
 * to the ccd.
 * This can be removed/altered if implementation changes.
 *
 */

void CrossValidationDriver::changeParameter(CyclicCoordinateDescent &ccd, int varianceIndex, double varianceValue) {
	if (varianceIndex == 0) {
		ccd.sigma2Beta = varianceValue;

	}
	if (varianceIndex == 1) {
		ccd.classHierarchyVariance = varianceValue;
	}
}

/* Author: tshaddox
 * This performs cross validation over one fold.  This is a way of
 * decomposing the actions of the greedy cross validation.
 */

double CrossValidationDriver::oneFoldCrossValidation(CyclicCoordinateDescent& ccd,
		AbstractSelector& selector,
		const CCDArguments& arguments,
		int i, vector<bsccs::real> weights,
		int step, int point) {

	int fold = i % arguments.fold;
	if (fold == 0) {
		selector.permute(); // Permute every full cross-validation rep
	}


	// Get this fold and update
	selector.getWeights(fold, weights);
	ccd.setWeights(&weights[0]);
	ccd.update(arguments.maxIterations, arguments.convergenceType, arguments.tolerance);  //tshaddox temporary comment out

	// Compute predictive loglikelihood for this fold
	selector.getComplement(weights);
	double logLikelihood = ccd.getPredictiveLogLikelihood(&weights[0]);
	return(logLikelihood);

}

/* Pseudocode for the algorithm:
	 * While the predicted log likelihood is not changing much...
	 * 	For each element sigma2Beta_i in the vector of variances:
	 * 		Fix all other elements in the vector of variances
	 *		Perform cross validation along the grid points for sigma2Beta_i
	 *		Select the value of sigma2Beta_i that produces maximal predicted log likelihood
	 *	If the number of iterations maximal is exceeded, stop.
	 */

void CrossValidationDriver::greedyDrive(CyclicCoordinateDescent& ccd,
		AbstractSelector& selector,
		const CCDArguments& arguments) {

	double storedMaxValue;
	bool finished = false;
	int count = 0;

	std::vector<double> varianceVector;
	int nVariances = 2;

	//initialize variances
	for (int i = 0; i < nVariances; i ++){
		varianceVector.push_back(computeGridPoint(0));
		changeParameter(ccd, i, computeGridPoint(0));
	}

	while (!finished) {

		if (count > 10) {
			finished = true;
		}
		// For each variance...
		for (int varianceIndex = 0; varianceIndex < varianceVector.size(); varianceIndex++) {

			std::vector<bsccs::real> weights;
			std::vector<double> points;
			std::vector<double> values;

			//Walk along the parameter to find the optimum
			for (int step = 0; step < gridSize; step++) {
				std::vector<double> predLogLikelihood;
				double point = computeGridPoint(step);
				varianceVector[varianceIndex] = point;
				changeParameter(ccd, varianceIndex, point); //set hyperprior to "point" value
				//Perform k fold cross validation
				for (int i = 0; i < arguments.foldToCompute; i++) {
					double logLikelihood = oneFoldCrossValidation(ccd, selector, arguments, i, weights, step, point);
					predLogLikelihood.push_back(logLikelihood);
				}

				double value = computePointEstimate(predLogLikelihood) /
						(double(arguments.foldToCompute) / double(arguments.fold));
				points.push_back(point);
				values.push_back(value);
			}
			// Report results
			double maxPoint = points[0];
			double maxValue = values[0];
			for(int k = 0; k < points.size(); k++) {
				if(values[k] > maxValue){
					maxPoint = points[k];
					maxValue = values[k];
				}
			}
			cout << endl;
			cout << endl;
			cout << "Max value is " << maxValue << endl;
			cout << "Max point is " << maxPoint << endl;
			changeParameter(ccd, varianceIndex, maxPoint);
			cout << endl;
			cout << endl;
			double testValue = storedMaxValue - maxValue;
			if (testValue*testValue < 10) {
				finished = true;
			}
			if (count + varianceIndex == 0 || maxValue > storedMaxValue){
				storedMaxValue = maxValue;
				cout << "---------------------  New Max = " << storedMaxValue << "-------------------" << endl;
			}
		}
		count++;
	}
}


void CrossValidationDriver::hierarchyDrive(CyclicCoordinateDescent& ccd,
		AbstractSelector& selector,
		const CCDArguments& arguments) {

	std::vector<bsccs::real> weights;
	std::vector<double> outerPoints;
	std::vector<double> innerPoints;
	std::vector<double> outerValues;
	std::vector<double> minValues;
	for (int outerStep = 0; outerStep < gridSize; outerStep++){
		std::vector<double> predLogLikelihoodOuter;
		double outerPoint = computeGridPoint(outerStep);
		ccd.classHierarchyVariance = outerPoint;
		cout << outerStep << " out of " << gridSize << endl;

		for (int step = 0; step < gridSize; step++) {

			std::vector<double> predLogLikelihood;
			double point = computeGridPoint(step);
			ccd.setHyperprior(point);

			for (int i = 0; i < arguments.foldToCompute; i++) {

				int fold = i % arguments.fold;
				if (fold == 0) {
					selector.permute(); // Permute every full cross-validation rep
				}

				// Get this fold and update
				selector.getWeights(fold, weights);
				ccd.setWeights(&weights[0]);
			//	std::cout << "Running at " << ccd.getPriorInfo() << " ";
				ccd.update(arguments.maxIterations, arguments.convergenceType, arguments.tolerance);  //tshaddox temporary comment out
				// Compute predictive loglikelihood for this fold
				selector.getComplement(weights);
				double logLikelihood = ccd.getPredictiveLogLikelihood(&weights[0]);
/*
				std::cout << "Grid-point #" << (step + 1) << " at " << point;
				std::cout << "\tFold #" << (fold + 1)
				          << " Rep #" << (i / arguments.fold + 1) << " pred log like = "
				          << logLikelihood << std::endl;
*/
				// Store value
				predLogLikelihood.push_back(logLikelihood);
			}

			double value = computePointEstimate(predLogLikelihood) /
					(double(arguments.foldToCompute) / double(arguments.fold));
			gridPoint.push_back(point);
			gridValue.push_back(value);
		}

		// Report results
		double maxPoint;
		double maxValue;
		double minValue;
		findMax(&maxPoint, &maxValue, &minValue);
		minValues.push_back(minValue);

		innerPoints.push_back(maxPoint);
		outerPoints.push_back(outerPoint);
		outerValues.push_back(maxValue);



		std::cout << std::endl;
		std::cout << "Maximum predicted log likelihood (" << maxValue << ") found at:" << std::endl;
		std::cout << "\t" << maxPoint << " (variance)" << std::endl;

		if (!arguments.useNormalPrior) {
			double lambda = convertVarianceToHyperparameter(maxPoint);
			std::cout << "\t" << lambda << " (lambda)" << std::endl;
		}

		//std:cout << std::endl;

	}
	double outerMaxPoint = outerPoints[0];
	double innerMaxPoint = innerPoints[0];
	double outerMaxValue = outerValues[0];
	for (int i = 0; i < outerPoints.size(); i++) {
		if (outerValues[i] > outerMaxValue) {
			outerMaxValue = outerValues[i];
			outerMaxPoint = outerPoints[i];
			innerMaxPoint = innerPoints[i];
		}
		cout << minValues[i] << endl;
	}
	std::cout << std::endl;
	std::cout << "Maximum predicted log likelihood (" << outerMaxValue << ") found at:" << std::endl;
	std::cout << "\t" << innerMaxPoint << " (drug variance) and at " << outerMaxPoint << " (class variance)" << std::endl;

}

void CrossValidationDriver::drive(
		CyclicCoordinateDescent& ccd,
		AbstractSelector& selector,
		const CCDArguments& arguments) {

	// TODO Check that selector is type of CrossValidationSelector

	std::vector<bsccs::real> weights;

	for (int step = 0; step < gridSize; step++) {

		std::vector<double> predLogLikelihood;
		double point = computeGridPoint(step);
		ccd.setHyperprior(point);

		for (int i = 0; i < arguments.foldToCompute; i++) {

			int fold = i % arguments.fold;
			if (fold == 0) {
				selector.permute(); // Permute every full cross-validation rep
			}

			// Get this fold and update
			selector.getWeights(fold, weights);
			ccd.setWeights(&weights[0]);
			std::cout << "Running at " << ccd.getPriorInfo() << " ";
			ccd.update(arguments.maxIterations, arguments.convergenceType, arguments.tolerance);  //tshaddox temporary comment out

			// Compute predictive loglikelihood for this fold
			selector.getComplement(weights);
			double logLikelihood = ccd.getPredictiveLogLikelihood(&weights[0]);
/*
			std::cout << "Grid-point #" << (step + 1) << " at " << point;
			std::cout << "\tFold #" << (fold + 1)
			          << " Rep #" << (i / arguments.fold + 1) << " pred log like = "
			          << logLikelihood << std::endl;
*/
			// Store value
			predLogLikelihood.push_back(logLikelihood);
		}

		double value = computePointEstimate(predLogLikelihood) /
				(double(arguments.foldToCompute) / double(arguments.fold));
		gridPoint.push_back(point);
		gridValue.push_back(value);
	}

	// Report results
	double maxPoint;
	double maxValue;
	double minValue;
	findMax(&maxPoint, &maxValue, &minValue);

	std::cout << std::endl;
	std::cout << "Maximum predicted log likelihood (" << maxValue << ") found at:" << std::endl;
	std::cout << "\t" << maxPoint << " (variance)" << std::endl;

	if (!arguments.useNormalPrior) {
		double lambda = convertVarianceToHyperparameter(maxPoint);
		std::cout << "\t" << lambda << " (lambda)" << std::endl;
	}
	std:cout << std::endl;
}


void CrossValidationDriver::findMax(double* maxPoint, double* maxValue, double* minValue) {

	*maxPoint = gridPoint[0];
	*maxValue = gridValue[0];
	for (int i = 1; i < gridPoint.size(); i++) {
		//cout << "MaxValues i =" << gridValue[i] << endl;
		if (gridValue[i] > *maxValue) {
			*maxPoint = gridPoint[i];
			*maxValue = gridValue[i];
		}
		if (gridValue[i] < *minValue){
			*minValue = gridValue[i];
		}
	}
	cout << "Max = " << *maxValue << endl;
	cout << "Min = " << *minValue << endl;
}
}

