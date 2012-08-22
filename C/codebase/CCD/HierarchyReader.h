/*
 * HierarchyReader.h
 *
 *  Created on: Jun 7, 2011
 *      Author: tshaddox
 */

#ifndef HIERARCHYREADER_H_
#define HIERARCHYREADER_H_

#include <iostream>
#include <fstream>

#include <vector>
#include <map>

using namespace std;

#include "CyclicCoordinateDescent.h"

//#define USE_DRUG_STRING

namespace bsccs {

#ifdef USE_DRUG_STRING
	typedef string DrugIdType; // TODO String do not get sorted in numerical order
#else
	typedef int DrugIdType;
#endif

class HierarchyReader: public CompressedIndicatorMatrix {
public:
	HierarchyReader();

	HierarchyReader(const char* fileName, map<int, int> &drugIdToIndex, CyclicCoordinateDescent **ccd);


private:


};

}

#endif /* HIERARCHYREADER_H_ */
