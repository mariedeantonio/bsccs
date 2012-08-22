/*
 * ccd.h
 *
 *  Created on: Sep 10, 2010
 *      Author: msuchard
 */

#ifndef CCD_H_
#define CCD_H_

#include <time.h>
#include <sys/time.h>

#include "CyclicCoordinateDescent.h"

using namespace bsccs;

struct CCDArguments {

	// Needed for fitting
	std::string inFileName;
	std::string outFileName;
	std::string hierarchyFileName; //tshaddox
	double classHierarchyVariance; //tshaddox
	double sigma2Beta; //tshaddox
	bool useGPU;
	bool useBetterGPU;
	int deviceNumber;
	double tolerance;
	double hyperprior;
	bool useNormalPrior;
	bool hyperPriorSet;
	int maxIterations;
	int convergenceType;
	long seed;

	// Needed for cross-validation
	bool doCrossValidation;
	double lowerLimit;
	double upperLimit;
	int fold;
	int foldToCompute;
	int gridSteps;
	std::string cvFileName;
	bool doFitAtOptimal;

	// Needed for boot-strapping
	bool doBootstrap;
	bool reportRawEstimates;
	bool BCaBootstrapCI;
	int replicates;
	std::string bsFileName;
};

void parseCommandLine(
		int argc,
		char* argv[],
		CCDArguments &arguments);

double initializeModel(
		InputReader** reader,
		CyclicCoordinateDescent** ccd,
		CCDArguments &arguments, map<int, int> &drugIdToIndex);

double fitModel(
		CyclicCoordinateDescent *ccd,
		CCDArguments &arguments);

double runCrossValidation(
		CyclicCoordinateDescent *ccd,
		InputReader *reader,
		CCDArguments &arguments);

double runBoostrap(
		CyclicCoordinateDescent *ccd,
		InputReader *reader,
		CCDArguments &arguments);

double calculateSeconds(
		const struct timeval &time1,
		const struct timeval &time2);

#endif /* CCD_H_ */



