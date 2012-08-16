/*
 * ImputeInputReader.cpp
 *
 *  Created on: Aug 13, 2012
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

#include "BBRImputeInputReader.h"

#define MAX_ENTRIES		1000000000

#define FORMAT_MATCH_1	"CONDITION_CONCEPT_ID"
#define MATCH_LENGTH_1	20
#define FORMAT_MATCH_2	"PID"
#define MATCH_LENGTH_2	3
#define MISSING_STRING_1	"NA"
#define MISSING_STRING_2	"NULL"
#define MISSING_LENGTH	-1
#define DELIMITER		" "

using namespace std;

BBRImputeInputReader::BBRImputeInputReader() : ImputeInputReader() { }

BBRImputeInputReader::~BBRImputeInputReader() { }

/**
 * Reads in a dense BBR data file with format:
 * Stratum,Outcome,X1 ...
 *
 * Assumes that file is sorted by 'Stratum'
 */

void BBRImputeInputReader::readFile(const char* fileName) {

	ifstream in(fileName);
	if (!in) {
		cerr << "Unable to open " << fileName << endl;
		exit(-1);
	}

	string line;
//	getline(in, line); // Read header and ignore

	int numCases = 0;
	vector<string> strVector;
	string outerDelimiter(DELIMITER);

	// Allocate a column for intercept
	real_vector* thisColumn = new real_vector();
	push_back(NULL, thisColumn, INDICATOR);
	int_vector* nullVector1 = new int_vector();
	int_vector* nullVector2 = new int_vector();
	entriesPresent.push_back(nullVector1);
	entriesAbsent.push_back(nullVector2);
	columnType.push_back("lr");
	nMissingPerColumn.push_back(0);
	drugMap.insert(make_pair(0,0));
	indexToDrugIdMap.insert(make_pair(0,0));

	int currentRow = 0;
	int maxCol = 0;
	while (getline(in, line) && (currentRow < MAX_ENTRIES)) {
		if (!line.empty()) {

			strVector.clear();
			split(strVector, line, outerDelimiter);

			nevents.push_back(1);
			numCases++;
			pid.push_back(numCases - 1);

			vector<string> thisCovariate;
			split(thisCovariate, strVector[0], ":");

			// Parse outcome entry
			real thisY = static_cast<real>(atof(thisCovariate[0].c_str()));
			y.push_back(thisY);

			// Fix offs for CLR
			offs.push_back(1);

			//Fill intercept
			data[0]->push_back(1.0);
			entriesPresent[0]->push_back(currentRow);

			// Parse covariates
			for (int i = 0; i < (int)strVector.size() - 1; ++i){
				thisCovariate.clear();
				split(thisCovariate, strVector[i + 1], ":");
				if((int)thisCovariate.size() == 2){
					DrugIdType drug = static_cast<DrugIdType>(atof(thisCovariate[0].c_str()));
					if(drugMap.count(drug) == 0){
						maxCol++;
						drugMap.insert(make_pair(drug,maxCol));
						indexToDrugIdMap.insert(make_pair(maxCol,drug));

						real_vector* thisColumn = new real_vector();
						push_back(NULL, thisColumn, INDICATOR);
						int_vector* nullVector1 = new int_vector();
						int_vector* nullVector2 = new int_vector();
						entriesPresent.push_back(nullVector1);
						entriesAbsent.push_back(nullVector2);
						columnType.push_back("lr");
						nMissingPerColumn.push_back(0);
					}
					int col = drugMap[drug];
					for(int j = (int)data[col]->size(); j < currentRow; j++){
						data[col]->push_back(0.0);
						entriesPresent[col]->push_back(j);
					}
					if(thisCovariate[1] == MISSING_STRING_1 || thisCovariate[1] == MISSING_STRING_2){
						data[col]->push_back(0.0);
						entriesAbsent[col]->push_back(currentRow);
						nMissingPerColumn[col]++;
					}
					else{
						real value = static_cast<real>(atof(thisCovariate[1].c_str()));
						data[col]->push_back(value);
						entriesPresent[col]->push_back(currentRow);
						if(value != 1.0 && value != 0.0)
						{
							columnType[col] = "ls";
							formatType[col] = DENSE;
						}
					}
				}
/*
				for(int j = col; j <= drug; j++){
					if(j > maxCol){
						real_vector* thisColumn = new real_vector();
						push_back(NULL, thisColumn, INDICATOR);
						int_vector* nullVector1 = new int_vector();
						int_vector* nullVector2 = new int_vector();
						entriesPresent.push_back(nullVector1);
						entriesAbsent.push_back(nullVector2);
						columnType.push_back("lr");
						nMissingPerColumn.push_back(0);
						maxCol++;
					}
					if(j != drug){
						data[j]->push_back(0.0);
						entriesPresent[j]->push_back(currentRow);
					}
					else{
						if(thisCovariate[1] == MISSING_STRING_1 || thisCovariate[1] == MISSING_STRING_2){
							data[j]->push_back(0.0);
							entriesAbsent[j]->push_back(currentRow);
							nMissingPerColumn[j]++;
						}
						else{
							real value = static_cast<real>(atof(thisCovariate[1].c_str()));
							data[j]->push_back(value);
							entriesPresent[j]->push_back(currentRow);
							if(value != 1.0 && value != 0.0)
							{
								columnType[j] = "ls";
								formatType[j] = DENSE;
							}
						}
					}
				}
				col = drug + 1;
				*/
			}
			currentRow++;
		}
	}

	nPatients = numCases;
	nCols = columns.size();
	nRows = currentRow;
	conditionId = "0";

	nCols_ = nCols;
	nRows_ = nRows;
	y_ = y;

	nevents.push_back(1); // Save last patient
	
	for(int j = 0; j < (int)columns.size(); j++){
		for(int i = (int)data[j]->size(); i < nRows; i++){
			data[j]->push_back(0.0);
			entriesPresent[j]->push_back(i);
		}
	}

	for(int j = 0; j < nCols; j++){
		if(formatType[j] == INDICATOR){
			if(columns[j])
				columns[j]->clear();
			columns[j] = new int_vector();
			for(int i = 0; i < nRows; i++)
			{
				if((int)data[j]->at(i) == 1)
					columns[j]->push_back(i);
			}
			data[j]->clear();
		}
	}

	int index = columns.size();

#ifndef MY_RCPP_FLAG
	cout << "ImputeInputReader" << endl;
	cout << "Read " << currentRow << " data lines from " << fileName << endl;
	cout << "Number of stratum: " << numCases << endl;
	cout << "Number of covariates: " << nCols << endl;
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
