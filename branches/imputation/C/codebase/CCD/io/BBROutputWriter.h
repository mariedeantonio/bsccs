/*
* BBROutputWriter.h
*
*  Created on: Aug 27, 2012
*      Author: Sushil Mittal
*/

#ifndef BBROUTPUTWRITER_H_
#define BBROUTPUTWRITER_H_

#include "ModelData.h"
#include "ImputationPolicy.h"

using namespace std;

class BBROutputWriter{
public:
	BBROutputWriter() {} ;
	virtual ~BBROutputWriter() {}

	void BBROutputWriter::writeFile(const char* fileName, ModelData* modelData, vector<int> reverseColIndices = vector<int>()) {
		ofstream out;
		out.open(fileName,ios::out);

		int nRows = modelData->getNumberOfRows();
		int nCols = modelData->getNumberOfColumns();

		if((int)reverseColIndices.size() == 0){
			for(int i = 0; i < nCols; i++)
				reverseColIndices.push_back(i);
		}

		map<DrugIdType,int> drugMap = modelData->getDrugMap();
		vector<real> y = modelData->getYVectorRef();
		vector<real> z = modelData->getZVectorRef();
		for(int i = 0; i < nRows; i++){
			vector<real> x(nCols,0.0);
			modelData->getDataRow(i,&x[0]);
			out << y[i];
			if((int)z.size() > 0)
				out << ":" << z[i];

			map<DrugIdType,int>::iterator it = drugMap.begin();
			for(int j = 0; j < nCols; j++, *it++){
				int drug = (*it).first;
				if(drug == 0)
					continue;
				int col = reverseColIndices[drugMap[drug]];
				if(x[col] != 0.0)
					out << " " << drug << ":" << x[col];
			}
			out << endl;
		}
	}
};

#endif /* BBROUTPUTWRITER_H_ */
