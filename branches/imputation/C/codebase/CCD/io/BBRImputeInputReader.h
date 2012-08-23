/*
 * ImputateInputReader.h
 *
 *  Created on: Aug 13, 2012
 *      Author: Sushil Mittal
 */

#ifndef BBRIMPUTEINPUTREADER_H_
#define BBRIMPUTEINPUTREADER_H_

#include "InputReader.h"
#include "ImputeInputReader.h"

class BBRImputeInputReader : public InputReader, public ImputeInputReader {
public:
	BBRImputeInputReader();
	virtual ~BBRImputeInputReader();

	virtual void readFile(const char* fileName);
//	void writeFile(const char* fileName);
};

#endif /* BBRIMPUTEINPUTREADER_H_ */
