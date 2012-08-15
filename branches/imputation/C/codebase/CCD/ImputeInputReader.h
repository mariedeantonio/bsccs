/*
 * ImputateInputReader.h
 *
 *  Created on: Jul 28, 2012
 *      Author: Sushil Mittal
 */

#ifndef IMPUTEINPUTREADER_H_
#define IMPUTEINPUTREADER_H_

#include "InputReader.h"

class Compare{
	vector<int>& _vec;
public:
	Compare(vector<int>& vec) : _vec(vec) {}
	bool operator()(size_t i, size_t j){
		return _vec[i] < _vec[j];
	}
};

class ImputeInputReader : public InputReader {
public:
	ImputeInputReader();
	virtual ~ImputeInputReader();

	virtual void readFile(const char* fileName) = 0;
	void sortColumns();
	vector<string> getColumnTypesToImpute();
	vector<int> getnMissingPerColumn();
	void setupDataForImputation(int col,vector<real>& weights);
	template <class T> void reindexVector(vector<T>& vec, vector<int> ind) {
		int n = (int) vec.size();
		vector<T> temp = vec;
		for(int i = 0; i < n; i++){
			vec[i] = temp[ind[i]];
		}
	}
	void resetData();
	void getSampleMeanVariance(int col, real* Xmean, real* Xvar);
	real* getDataRow(int row, real* x);
	void updateColumnVector(int col, vector<int> appendY);
protected:
	vector<string> columnType;
	vector<int> nMissingPerColumn;
	vector<int> colIndices;
	vector<real> y_;
	vector<int_vector*> entriesAbsent;
	vector<int_vector*> entriesPresent;

	int nCols_;
	int nRows_;
};

#endif /* IMPUTEINPUTREADER_H_ */
