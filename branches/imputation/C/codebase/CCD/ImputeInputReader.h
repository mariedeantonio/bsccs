/*
 * ImputateInputReader.h
 *
 *  Created on: Jul 28, 2012
 *      Author: Sushil Mittal
 */

#ifndef IMPUTEINPUTREADER_H_
#define IMPUTEINPUTREADER_H_

#include "InputReader.h"

class ImputeInputReader : public InputReader {
public:
	ImputeInputReader();
	virtual ~ImputeInputReader();

	virtual void readFile(const char* fileName);
	vector<string> getColumnTypesToImpute();
	vector<int> getnMissingPerColumn();
	void setupData(int col,vector<real>& weights);
	template <class T>
	void reindexVector(vector<T>& vec, vector<int> ind);
private:
	vector<string> columnType;
	vector<int> nMissingPerColumn;
	vector<int> colIndices;
	vector<real_vector*> data_;
	vector<int_vector*> entriesAbsent;
	vector<int_vector*> entriesPresent;

	int nPatients_;
	int nCols_;
	int nRows_;
};

#endif /* IMPUTEINPUTREADER_H_ */
