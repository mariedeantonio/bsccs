/*
 * ImputateInputReader.h
 *
 *  Created on: Aug 13, 2012
 *      Author: Sushil Mittal
 */

#ifndef BBRIMPUTEINPUTREADER_H_
#define BBRIMPUTEINPUTREADER_H_

#include "ImputeInputReader.h"

class BBRImputeInputReader : public ImputeInputReader {
public:
	BBRImputeInputReader();
	virtual ~BBRImputeInputReader();

	virtual void readFile(const char* fileName);
};

#endif /* BBRIMPUTEINPUTREADER_H_ */
