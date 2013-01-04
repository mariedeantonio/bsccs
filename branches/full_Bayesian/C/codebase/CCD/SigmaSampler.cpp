/*
 * SigmaCalculate.cpp
 *
 *  Created on: Aug 9, 2012
 *      Author: trevorshaddox
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <map>
#include <time.h>
#include <set>

#include <SigmaSampler.h>

#include <boost/random.hpp>
#include <boost/random/gamma_distribution.hpp>

using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;


namespace bsccs {


	SigmaSampler::SigmaSampler(){

	}

	SigmaSampler::~SigmaSampler(){

	}

	void SigmaSampler::sampleSigma(Parameter * SigmaSquared, Parameter * BetaValues, boost::mt19937& rng){
		//TODO fix passing in rng as boost::mt19937

		// tau | BetaVector ~ gamma(alpha + N/2, Beta + (1/2)(SUM(beta_i - mu)^2)
		// prior: tau ~ gamma(alpha, beta)

		double SigmaParameter_alpha = 2;
		double SigmaParameter_beta = 4;

		double BetaMinusMu = 0;
		double Mu = 0;

		for (int i = 0; i < BetaValues->getSize(); i++) {
			Mu += BetaValues->get(i) / BetaValues->getSize();
		}

		for (int j = 0; j < BetaValues->getSize(); j++ ) {
			BetaMinusMu += (BetaValues->get(j) - Mu)*(BetaValues->get(j) - Mu); //
		}

		const double shape = SigmaParameter_alpha + BetaValues->getSize() / 2;
		double scale = SigmaParameter_beta + BetaMinusMu / 2;

		boost::gamma_distribution<> gd( shape );
		boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma( rng, gd );

		double newValue = 1/scale*var_gamma();

		SigmaSquared->set(0, newValue);
		//SigmaSquared->logParameter();

	}



}

