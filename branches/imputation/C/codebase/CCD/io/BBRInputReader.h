/*
* BBRInputReader.h
*
*  Created on: Aug 13, 2012
*      Author: Sushil Mittal
*/

#ifndef BBRINPUTREADER_H_
#define BBRINPUTREADER_H_

#include "InputReader.h"
#include "ImputationPolicy.h"

using namespace std;

#define MAX_ENTRIES		1000000000

#define FORMAT_MATCH_1	"CONDITION_CONCEPT_ID"
#define MATCH_LENGTH_1	20
#define FORMAT_MATCH_2	"PID"
#define MATCH_LENGTH_2	3
#define MISSING_STRING_1	"NA"
#define MISSING_STRING_2	"NULL"
#define MISSING_LENGTH	-1
#define DELIMITER		" "

template <class ImputationPolicy>
class BBRInputReader : public InputReader{
public:
	BBRInputReader():InputReader() {
		imputePolicy = new ImputationPolicy();
	}
	virtual ~BBRInputReader() {}

	virtual void readFile(const char* fileName){

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
		modelData->push_back(NULL, thisColumn, DENSE);
		int_vector* nullVector = new int_vector();
		imputePolicy->push_back(nullVector,0);
		modelData->drugMap.insert(make_pair(0,0));
		modelData->indexToDrugIdMap.insert(make_pair(0,0));

		int currentRow = 0;
		modelData->nRows = 0;
		int maxCol = 0;
		while (getline(in, line) && (currentRow < MAX_ENTRIES)) {
			if (!line.empty()) {

				strVector.clear();
				split(strVector, line, outerDelimiter);

				modelData->nevents.push_back(1);
				numCases++;
				modelData->pid.push_back(numCases - 1);

				vector<string> thisCovariate;
				split(thisCovariate, strVector[0], ":");

				// Parse outcome entry
				real thisY = static_cast<real>(atof(thisCovariate[0].c_str()));
				modelData->y.push_back(thisY);

				// Parse censoring index entry if any
				if(thisCovariate.size() == 2){
					real thisZ = static_cast<real>(atof(thisCovariate[1].c_str()));
					modelData->z.push_back(thisZ);
				}

				// Fix offs for CLR
				modelData->offs.push_back(1);

				//Fill intercept
				modelData->data[0]->push_back(1.0);

				// Parse covariates
				for (int i = 0; i < (int)strVector.size() - 1; ++i){
					thisCovariate.clear();
					split(thisCovariate, strVector[i + 1], ":");
					if((int)thisCovariate.size() == 2){
						DrugIdType drug = static_cast<DrugIdType>(atof(thisCovariate[0].c_str()));
						if(modelData->drugMap.count(drug) == 0){
							maxCol++;
							modelData->drugMap.insert(make_pair(drug,maxCol));
							modelData->indexToDrugIdMap.insert(make_pair(maxCol,drug));

							int_vector* colInds = new int_vector();
							real_vector* thisColumn = new real_vector();
							modelData->push_back(colInds, thisColumn, INDICATOR);
							int_vector* nullVector = new int_vector();
							imputePolicy->push_back(nullVector,0);
						}
						int col = modelData->drugMap[drug];
						if(modelData->formatType[col] == DENSE){
							for(int j = (int)modelData->data[col]->size(); j < currentRow; j++){
								modelData->data[col]->push_back(0.0);
							}
						}
						if(thisCovariate[1] == MISSING_STRING_1 || thisCovariate[1] == MISSING_STRING_2){
							if(modelData->formatType[col] == DENSE){
								modelData->data[col]->push_back(0.0);
							}
							imputePolicy->push_back(col,currentRow);
						}
						else{
							real value = static_cast<real>(atof(thisCovariate[1].c_str()));
							if(modelData->formatType[col] == DENSE){
								modelData->data[col]->push_back(value);
							}
							else if(modelData->formatType[col] == INDICATOR){
								if(value != 1.0 && value != 0.0){
									if(modelData->columns[col]->size() > 0){
										modelData->convertColumnToDense(col);
									}
									else{
										delete modelData->columns[col];
										modelData->columns[col] = NULL;
									}
									for(int j = (int)modelData->data[col]->size(); j < currentRow; j++){
										modelData->data[col]->push_back(0.0);
									}
									modelData->data[col]->push_back(value);
									modelData->formatType[col] = DENSE;
								}
								else{
									if(value == 1.0){
										modelData->columns[col]->push_back(currentRow);
									}
								}
							}
						}
					}
				}
				currentRow++;
				modelData->nRows++;
			}
		}

		modelData->nPatients = numCases;
		modelData->nCols = modelData->columns.size();
		modelData->conditionId = "0";

		imputePolicy->setParams(modelData->y, modelData->nCols);

		modelData->nevents.push_back(1); // Save last patient

		for(int j = 0; j < modelData->nCols; j++){
			if(modelData->formatType[j] == DENSE){
				for(int i = (int)modelData->data[j]->size(); i < modelData->nRows; i++){
					modelData->data[j]->push_back(0.0);
				}
			}
		}
		int index = modelData->columns.size();

#ifndef MY_RCPP_FLAG
		cout << "ImputeInputReader" << endl;
		cout << "Read " << currentRow << " data lines from " << fileName << endl;
		cout << "Number of stratum: " << numCases << endl;
		cout << "Number of covariates: " << modelData->nCols << endl;
#endif

	}

	ImputationPolicy* getImputationPolicy(){
		return imputePolicy;
	}

protected:
	ImputationPolicy* imputePolicy;
};

#endif /* BBRINPUTREADER_H_ */
