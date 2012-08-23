/*
 * ImputateInputReader.h
 *
 *  Created on: Jul 28, 2012
 *      Author: Sushil Mittal
 */

#ifndef IMPUTEINPUTREADER_H_
#define IMPUTEINPUTREADER_H_

#include "io/InputReader.h"

class Compare{
	vector<int>& _vec;
public:
	Compare(vector<int>& vec) : _vec(vec) {}
	bool operator()(size_t i, size_t j){
		return _vec[i] < _vec[j];
	}
};

class ImputeInputReader{
public:
	ImputeInputReader();
	virtual ~ImputeInputReader();

	void sortColumns(vector<real_vector*>& data, vector<int_vector*>& columns, vector<FormatType>& formatType, int& nCols);
	vector<string> getColumnTypesToImpute();
	vector<int> getnMissingPerColumn();
	void setupDataForImputation(int col, vector<real>& weights, vector<real>& y, vector<real_vector*>& data, vector<int_vector*>& columns, vector<FormatType>& formatType, int& nRows, int& nCols);
	template <class T> void reindexVector(vector<T>& vec, vector<int> ind) {
		int n = (int) vec.size();
		vector<T> temp = vec;
		for(int i = 0; i < n; i++){
			vec[i] = temp[ind[i]];
		}
	}
	void resetParams(vector<real> y, int& nCols);
	void resetData(vector<int_vector*>& columns, vector<FormatType>& formatType, int& nCols);
	void getSampleMeanVariance(int col, real* Xmean, real* Xvar, vector<real_vector*>& data, vector<int_vector*>& columns, vector<FormatType>& formatType, int& nRows);
	real* getDataRow(int row, real* x, vector<real_vector*>& data, vector<int_vector*>& columns, vector<FormatType>& formatType, int& nCols);
	void updateColumnVector(int col, vector<int> appendY, vector<int_vector*>& columns);
protected:
	vector<string> columnType;
	vector<int> nMissingPerColumn;
	vector<int> colIndices;
	vector<int> reverseColIndices;
	vector<real> y_;
	vector<int_vector*> entriesAbsent;

	int nCols_;
};

#endif /* IMPUTEINPUTREADER_H_ */
