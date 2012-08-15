/*
 * ImputateInputReader.h
 *
 *  Created on: Aug 13, 2012
 *      Author: Sushil Mittal
 */

#ifndef CSVIMPUTEINPUTREADER_H_
#define CSVIMPUTEINPUTREADER_H_

#include "ImputeInputReader.h"

class CSVImputeInputReader : public ImputeInputReader {
public:
	CSVImputeInputReader();
	virtual ~CSVImputeInputReader();

	virtual void readFile(const char* fileName);
};

#endif /* CSVIMPUTEINPUTREADER_H_ */
