/*
 * BootstrapDriver.h
 *
 *  Created on: Nov 17, 2010
 *      Author: msuchard
 */

#ifndef BOOTSTRAPDRIVER_H_
#define BOOTSTRAPDRIVER_H_

#include <vector>

#include "AbstractDriver.h"
#include "InputReader.h"

typedef std::vector<bsccs::real> rvector;
typedef std::vector<rvector*> rarray;
typedef	rarray::iterator rarrayIterator;

namespace bsccs {

class BootstrapDriver : public AbstractDriver {
public:
	BootstrapDriver(
			int inReplicates,
			InputReader* inReader);

	virtual ~BootstrapDriver();

	virtual void drive(
			CyclicCoordinateDescent& ccd,
			AbstractSelector& selector,
			const CCDArguments& arguments);

	virtual void logResults(const CCDArguments& arguments);

private:
	double numeratorForBCa; // =sum(thetaHat - theta_i)^3 for all bootstrap samples
	double denominatorForBCa; // = 6*(sum(thetaHat - theta_i)^2)^3/2 for all bootstrap samples
	double BCaQ;
	const int replicates;
	InputReader* reader;
	const int J;
	rarray estimates;
};

}
#endif /* BOOTSTRAPDRIVER_H_ */
