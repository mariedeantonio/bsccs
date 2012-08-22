/*
 * BootstrapDriver.cpp
 *
 *  Created on: Nov 17, 2010
 *      Author: msuchard
 */

#include <iostream>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <cstdlib>
#include <math.h>

#include <boost/math/distributions/normal.hpp>

#include "BootstrapDriver.h"
#include "AbstractSelector.h"

using boost::math::normal;
namespace bsccs {

BootstrapDriver::BootstrapDriver(
		int inReplicates,
		InputReader* inReader) : replicates(inReplicates), reader(inReader),
		J(inReader->getNumberOfColumns()) {

	// Set-up storage for bootstrap estimates
	estimates.resize(J);
	int count = 0;
	for (rarrayIterator it = estimates.begin(); it != estimates.end(); ++it) {
		*it = new rvector();
	}
}

BootstrapDriver::~BootstrapDriver() {
	for (rarrayIterator it = estimates.begin(); it != estimates.end(); ++it) {
		if (*it) {
			delete *it;
		}
	}
}

void BootstrapDriver::drive(
		CyclicCoordinateDescent& ccd,
		AbstractSelector& selector,
		const CCDArguments& arguments) {

	// TODO Make sure that selector is type-of BootstrapSelector
	std::vector<bsccs::real> weights;
	struct timeval time1, time2;

	for (int j = 0; j < J; ++j) {
		estimates[j]->push_back(ccd.getBeta(j));
	}

	for (int step = 0; step < replicates; step++) {
		gettimeofday(&time1, NULL);
		selector.permute();
		selector.getWeights(0, weights);

		ccd.setWeights(&weights[0]);

		std::cout << std::endl << "Running replicate #" << (step + 1) << std::endl;
		// Run CCD using a warm start
		ccd.update(arguments.maxIterations, arguments.convergenceType, arguments.tolerance);
		gettimeofday(&time2, NULL);
		cout << "Time difference 1 = " << calculateSeconds(time1, time2) << endl;

		// Store point estimates
		for (int j = 0; j < J; ++j) {
			estimates[j]->push_back(ccd.getBeta(j));
		}
		gettimeofday(&time2, NULL);
		cout << "Time difference 2 = " << calculateSeconds(time1, time2) << endl;
	}
}

void BootstrapDriver::logResults(const CCDArguments& arguments) {

	ofstream outLog(arguments.bsFileName.c_str());
	if (!outLog) {
		cerr << "Unable to open log file: " << arguments.bsFileName << endl;
		exit(-1);
	}

	map<int, DrugIdType> drugMap = reader->getDrugNameMap();

	string sep(","); // TODO Make option

	if (!arguments.reportRawEstimates) {
	//	outLog << "drug" << sep << "mean" << sep << "var" << sep << "lower" << sep <<
	//		"upper" << sep << "prob0" << endl;
	}


	for (int j = 0; j < J; ++j) {

		outLog << drugMap[j] << sep;
		if (arguments.reportRawEstimates) {
			ostream_iterator<bsccs::real> output(outLog, sep.c_str());
			copy(estimates[j]->begin(), estimates[j]->end(), output);
			outLog << endl;
		} else {
			numeratorForBCa = 0;
			denominatorForBCa = 0;
			double nReplicates = 0;
			bsccs::real mean = 0.0;
			bsccs::real var = 0.0;
			bsccs::real prob0 = 0.0;
			double p = 0;
			bsccs::real thetaHat = *estimates[j]->begin();
			for (rvector::iterator it = estimates[j]->begin()+1; it != estimates[j]->end(); ++it) {
				numeratorForBCa += pow(thetaHat - *it,3);
				denominatorForBCa += pow(thetaHat - *it, 2);
				if (*it < thetaHat) {
					p++;
				}
				mean += *it;
				var += *it * *it;
				if (*it == 0.0) {
					prob0 += 1.0;
				}
				nReplicates ++;
			}
			normal s;
			//cout << "pre 3/2 denom = " << denominatorForBCa << endl;
			denominatorForBCa = 6*pow(denominatorForBCa, 3.0/2.0);

			//cout << "num = " << numeratorForBCa << endl;
			//cout << "denominator = " << denominatorForBCa << endl;

			double a = numeratorForBCa / denominatorForBCa;

			//cout << "a = " << a << endl;

			//cout << "p = " << p << endl;
			//cout << "nReplicates = " << nReplicates << endl;

			double zValue_lower = quantile(s,0.025);
			double zValue_upper = quantile(s,0.975);

			double bu = quantile(s, p / nReplicates);
			double bl = quantile(s, (nReplicates - p) / nReplicates);

			//cout << "bl = " << bl << endl;
			//cout << "bu = " << bu << endl;
			double Q_upper = cdf(s, bu - ((zValue_lower - bu) / (1+a*(zValue_lower-bu))));
			double Q_lower = cdf(s, bl - ((zValue_upper - bl) / (1+a*(zValue_upper-bl))));

			//cout << "Q_upper = " << Q_upper << endl;
			//cout << "Q_lower = " << Q_lower << endl;


			bsccs::real size = static_cast<bsccs::real>(estimates[j]->size());
			mean /= size;
			var = (var / size) - (mean * mean);
			prob0 /= size;

			sort(estimates[j]->begin(), estimates[j]->end());
			int offsetLower = static_cast<int>(size * 0.025);
			int offsetUpper = static_cast<int>(size * 0.975);

			int offsetLowerBCa = static_cast<int>(size * Q_lower);
			int offsetUpperBCa = static_cast<int>(size * Q_upper);

			// BCa Bootstrap Intervals

			bsccs::real BCaUpper = *(estimates[j]->begin() + offsetUpperBCa);
			cout << "BCaUpper = " << BCaUpper << endl;

			bsccs::real BCaLower = *(estimates[j]->begin() + offsetLowerBCa);
			cout << "BCaLower = " << BCaLower << endl;

			bsccs::real lower = *(estimates[j]->begin() + offsetLower);
			bsccs::real upper = *(estimates[j]->begin() + offsetUpper);

			//cout << "percentile Upper = " << upper << endl;
			//cout << "percentile Lower = " << lower << endl;


			if (arguments.BCaBootstrapCI) {
				outLog << mean << sep << var << sep << BCaLower << sep << BCaUpper << sep << prob0 << endl;
			} else {
				outLog << mean << sep << var << sep << lower << sep << upper << sep << prob0 << endl;
			}
		}

	}
	outLog.close();
}

}
