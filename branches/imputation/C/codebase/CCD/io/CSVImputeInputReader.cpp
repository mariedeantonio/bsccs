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

#include "CSVImputeInputReader.h"

#define MAX_ENTRIES		1000000000

#define FORMAT_MATCH_1	"CONDITION_CONCEPT_ID"
#define MATCH_LENGTH_1	20
#define FORMAT_MATCH_2	"PID"
#define MATCH_LENGTH_2	3
#define MISSING_STRING_1	"NA"
#define MISSING_STRING_2	"NULL"
#define MISSING_LENGTH	-1
#define DELIMITER		","

using namespace std;

CSVImputeInputReader::CSVImputeInputReader() : InputReader() { }

CSVImputeInputReader::~CSVImputeInputReader() { }

/**
 * Reads in a dense CSV data file with format:
 * Stratum,Outcome,X1 ...
 *
 * Assumes that file is sorted by 'Stratum'
 */


void CSVImputeInputReader::readFile(const char* fileName) {
	ifstream in(fileName);
	if (!in) {
		cerr << "Unable to open " << fileName << endl;
		exit(-1);
	}

	string line;
	getline(in, line); // Read header and ignore

	int numCases = 0;
	int numCovariates = MISSING_LENGTH;
	string currentStratum = MISSING_STRING_1;
	int numEvents = 0;

	vector<string> strVector;
	string outerDelimiter(DELIMITER);

	int currentRow = 0;
	while (getline(in, line) && (currentRow < MAX_ENTRIES)) {
		if (!line.empty()) {

			strVector.clear();
			split(strVector, line, outerDelimiter);

			// Make columns
			if (numCovariates == MISSING_LENGTH) {
				numCovariates = strVector.size() - 2;
				for (int i = 0; i < numCovariates; ++i) {
					real_vector* thisColumn = new real_vector();
					modelData->push_back(NULL, thisColumn, INDICATOR);
					int_vector* nullVector1 = new int_vector();
					int_vector* nullVector2 = new int_vector();
					entriesAbsent.push_back(nullVector2);
					columnType.push_back("lr");
				}
				nMissingPerColumn.resize(numCovariates,0);
			} else if (numCovariates != strVector.size() - 2) {
				cerr << "All rows must be the same length" << endl;
				exit(-1);
			}

			// Parse stratum (pid)
			string unmappedStratum = strVector[0];
			if (unmappedStratum != currentStratum) { // New stratum, ASSUMES these are sorted
				if (currentStratum != MISSING_STRING_1 && currentStratum != MISSING_STRING_2) { // Skip first switch
					modelData->nevents.push_back(1);
					numEvents = 0;
				}
				currentStratum = unmappedStratum;
				numCases++;
			}
			modelData->pid.push_back(numCases - 1);

			// Parse outcome entry
//			real thisY;
//			stringstream(strVector[1]) >> thisY;
			real thisY = static_cast<real>(atof(strVector[1].c_str()));
 			numEvents += thisY;
			modelData->y.push_back(thisY);

			// Fix offs for CLR
			modelData->offs.push_back(1);

			// Parse covariates
			for (int i = 0; i < numCovariates; ++i) {
//				real value;
//				istringstream(strVector[2 + i]) >> value;
				if(strVector[2 + i] == MISSING_STRING_1 || strVector[2 + i] == MISSING_STRING_2){
					modelData->data[i]->push_back(0);
					entriesAbsent[i]->push_back(currentRow);
						nMissingPerColumn[i]++;
				}
				else{
					real value = static_cast<real>(atof(strVector[2 + i].c_str()));
					modelData->data[i]->push_back(value);
					if(value != 1.0 && value != 0.0)
					{
						columnType[i] = "ls";
						modelData->formatType[i] = DENSE;
					}
				}
			}
			currentRow++;
		}
	}

	modelData->nPatients = numCases;
	modelData->nCols = modelData->columns.size();
	modelData->nRows = currentRow;
	modelData->conditionId = "0";

	nCols_ = modelData->nCols;
	y_ = modelData->y;

	modelData->nevents.push_back(1); // Save last patient

	for(int j = 0; j < modelData->nCols; j++){
		if(modelData->formatType[j] == INDICATOR){
			if(modelData->columns[j])
				modelData->columns[j]->clear();
			modelData->columns[j] = new int_vector();
			for(int i = 0; i < modelData->nRows; i++)
			{
				if((int)modelData->data[j]->at(i) == 1)
					modelData->columns[j]->push_back(i);
			}
			modelData->data[j]->clear();
		}
	}

	int index = modelData->columns.size();

#ifndef MY_RCPP_FLAG
	cout << "ImputeInputReader" << endl;
	cout << "Read " << currentRow << " data lines from " << fileName << endl;
	cout << "Number of stratum: " << numCases << endl;
	cout << "Number of covariates: " << numCovariates << endl;
#endif

#if 0
	for (int i = 0; i < nCols; ++i) {
		printColumn(i);
	}

	cerr << "PIDs ";
	printVector(&(pid[0]), static_cast<int>(pid.size()));

	cerr << "Ys ";
	printVector(eta.data(), eta.size());

	cerr << "nEvents ";
	printVector(nevents.data(), nevents.size());
#endif
}
/*
void CSVImputeInputReader::writeFile(const char* fileName){

}
*/