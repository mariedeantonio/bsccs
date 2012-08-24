/*
 * ImputateInputReader.h
 *
 *  Created on: Jul 28, 2012
 *      Author: Sushil Mittal
 */

#ifndef ImputationHelper_H_
#define ImputationHelper_H_

#include "io/InputReader.h"

class Compare{
	vector<int>& _vec;
public:
	Compare(vector<int>& vec) : _vec(vec) {}
	bool operator()(size_t i, size_t j){
		return _vec[i] < _vec[j];
	}
};

class ImputationHelper{
public:
	ImputationHelper();
	virtual ~ImputationHelper();

	void sortColumns();
	void push_back(int_vector* vecAbsent,int valMissing);
	void push_back(int col, int indAbsent);
	vector<int> getnMissingPerColumn();
	vector<int> getSortedColIndices();
	void setWeightsForImputation(int col, vector<real>& weights, int nRows);
	void setParams(vector<real> y, int nCols);
	int getNumberOfColumns();
	vector<real> getYVector();
	void getMissingEntries(int col, vector<int>& missing);
	void getSampleMeanVariance(int col, real* Xmean, real* Xvar, real* dataVec, int* columnVec, FormatType formatType, int nRows, int nEntries);
protected:
	vector<int> nMissingPerColumn;
	vector<int> colIndices;
	vector<int> reverseColIndices;
	vector<real> y_;
	vector<int_vector*> missingEntries;

	int nCols_;
};


class NoImputation{
public:
	NoImputation();
	virtual ~NoImputation();

	void sortColumns() {}
	void push_back(int_vector* vecAbsent,int valMissing) {}
	void push_back(int col, int indAbsent) {}
	vector<int> getnMissingPerColumn() {}
	vector<int> getSortedColIndices() {}
	void setWeightsForImputation(int col, vector<real>& weights, int nRows) {}
	void setParams(vector<real> y, int nCols) {}
	int getNumberOfColumns() {}
	vector<real> getYVector() {}
	void getMissingEntries(int col, vector<int>& missing) {}
	void getSampleMeanVariance(int col, real* Xmean, real* Xvar, real* dataVec, int* columnVec, FormatType formatType, int nRows, int nEntries) {}
};

#endif /* ImputationHelper_H_ */
