/*
 * ImputeInputReader.cpp
 *
 *  Created on: Jul 28, 2012
 *      Author: Sushil Mittal
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#include "ImputeInputReader.h"

#define MAX_ENTRIES		1000000000

#define FORMAT_MATCH_1	"CONDITION_CONCEPT_ID"
#define MATCH_LENGTH_1	20
#define FORMAT_MATCH_2	"PID"
#define MATCH_LENGTH_2	3
#define MISSING_STRING	"NA"
#define MISSING_LENGTH	-1
#define DELIMITER		","

using namespace std;

ImputeInputReader::ImputeInputReader() : InputReader() { }

ImputeInputReader::~ImputeInputReader() { }

/**
 * Reads in a dense CSV data file with format:
 * Stratum,Outcome,X1 ...
 *
 * Assumes that file is sorted by 'Stratum'
 */

template <class InputIterator1, class InputIterator2>
int set_intersection(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2)
{
	int result = 0;
	while (first1!=last1 && first2!=last2)
	{
		if(*first1 < *first2) 
			++first1;
		else if(*first2 < *first1) 
			++first2;
		else{ 
			result++;
			first1++; 
			first2++; 
		}
	}
	return result;
}

template <class InputIterator1, class InputIterator2>
int set_difference(InputIterator1 first1, InputIterator1 last1,	InputIterator2 first2, InputIterator2 last2)
{
	int result = last1 - first1;
	while (first1!=last1 && first2!=last2)
	{
		if(*first1 < *first2)
			++first1;
		else if(*first2 < *first1) 
			++first2;
		else{ 
			result--;
			first1++; 
			first2++; 
		}
	}
	return result;
}
void ImputeInputReader::sortColumns(){
	for(int i = 0; i < nCols; i++)
		colIndices.push_back(i);
	sort(colIndices.begin(),colIndices.end(),Compare(nMissingPerColumn));

	reverseColIndices.resize(nCols,0);
	for(int i = 0; i < nCols; i++)
		reverseColIndices[colIndices[i]] = i;

	reindexVector(nMissingPerColumn,colIndices);
	reindexVector(columnType,colIndices);
	reindexVector(formatType,colIndices);

	vector<real_vector*> tempData = data;
	vector<int_vector*> colData = columns;
	vector<int_vector*> entriesAbsent_ = entriesAbsent;
	for(int i = 0; i < nCols; i++)
	{
		data[i] = tempData[colIndices[i]];
		columns[i] = colData[colIndices[i]];
		entriesAbsent[i] = entriesAbsent_[colIndices[i]];
	}
}
vector<string> ImputeInputReader::getColumnTypesToImpute(){
	return columnType;
}

vector<int> ImputeInputReader::getnMissingPerColumn(){
	return nMissingPerColumn;
}

void ImputeInputReader::setupDataForImputation(int col, vector<real>& weights){
	y.clear();
	weights.clear();
	y.resize(nRows,0.0);
	weights.resize(nRows,1.0);
	if(formatType[col] == DENSE){
		for(int i = 0; i < nRows; i++)
			y[i] = data[col]->at(i);
	}
	else{
		for(int i = 0; i < (int)columns[col]->size(); i++)
			y[(int)columns[col]->at(i)] = 1.0;
	}
	for(int i = 0; i < (int)entriesAbsent[col]->size(); i++)
		weights[(entriesAbsent[col])->at(i)] = 0.0;

	nCols = col;

/*
	FILE* fp = fopen("Release/DATA.txt","w");
	for(int i = 0; i < nRows; i++)
	{
		for(int j = 0; j < nCols; j++)
		{
			fprintf(fp,"%f ",data[j]->at(i));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
*/
}

void ImputeInputReader::resetParams(){
	y.clear();
	nCols = nCols_;
	y = y_;
}

void ImputeInputReader::resetData(){
	for(int j = 0; j < nCols; j++){
		if(formatType[j] == INDICATOR){
			int lastit = 0;
			vector<int>::iterator it1 = entriesAbsent[j]->begin();
			vector<int>::iterator it2 = columns[j]->begin();
			while(it1 < entriesAbsent[j]->end() && it2 < columns[j]->end()){
				if(*it1 < *it2)
					it1++;
				else if(*it2 < *it1){
					it2++;
				}
				else{
					columns[j]->erase(it2);
					it2 = columns[j]->begin() + lastit;
				}
			}
		}
	}
}

void ImputeInputReader::getSampleMeanVariance(int col, real* Xmean, real* Xvar){
	for(int j = 0; j < col; j++){
		real sumx2 = 0.0;
		real sumx = 0.0;
		int n = nRows - (int)entriesAbsent[col]->size();
		if(formatType[j] == DENSE) {
			int ind = 0;
			for(int i = 0; i < nRows; i++){
				if(ind < (int)entriesAbsent[col]->size()){
					if(i == entriesAbsent[col]->at(ind)){
						ind++;
					}
				}
				real xi = data[j]->at(i);
				sumx2 += xi * xi;
				sumx += xi;
			}
		}
		else{
			sumx = set_difference(columns[j]->begin(), columns[j]->end(), entriesAbsent[col]->begin(), entriesAbsent[col]->end());
			sumx2 = sumx;
		}
		Xmean[j] = sumx/n;
		Xvar[j] =  (sumx2 - Xmean[j]*Xmean[j]*n)/(n-1);
	}
}

real* ImputeInputReader::getDataRow(int row, real* x){
	for(int j = 0; j < nCols; j++)
	{
		if(formatType[j] == DENSE)
			x[j] = data[j]->at(row);
		else{
			x[j] = 0.0;
			for(int i = 0; i < (int)columns[j]->size(); i++){
				if(columns[j]->at(i) == row){
					x[j] = 1.0;
					break;
				}
				else if(columns[j]->at(i) > row)
					break;
			}
		}

	}

	return x;
}

void ImputeInputReader::updateColumnVector(int col, vector<int> appendY){
	int lastit = 0;
	int p = columns[col]->size() + appendY.size();
	for(int i = 0; i < (int)appendY.size(); i++)
	{
		vector<int>::iterator it = columns[col]->begin() + lastit;
		if(columns[col]->size() > 0){
			while(*it < appendY[i]){
				it++;
				lastit++;
			}
		}
		columns[col]->insert(it,appendY[i]);
	}
}