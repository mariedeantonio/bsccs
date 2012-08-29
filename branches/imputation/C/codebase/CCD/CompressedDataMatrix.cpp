/*
 * CompressedDataMatrix.cpp
 *
 *  Created on: May-June, 2010
 *      Author: msuchard
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "CompressedDataMatrix.h"

CompressedDataMatrix::CompressedDataMatrix() {
	// Do nothing
}

//CompressedDataMatrix::CompressedDataMatrix(const char* fileName) {
//
//	ifstream in(fileName);
//	if (!in) {
//		cerr << "Unable to open " << fileName << endl;
//		exit(-1);
//	}
//
//	// Read header line
//	char buffer[256];
//	in.getline(buffer, 256);
//
//	// Read matrix dimensions
//	in >> nRows >> nCols >> nEntries;
//
////	// Allocate some memory
////	columns = std::vector<int_vector>(nCols);
////	for (int j = 0; j < nCols; j++) {
////		columns[j] = int_vector(); // Create empty list
////	}
//	allocateMemory(nCols);
//
//	// Read each matrix entry
//	for (int k = 0; k < nEntries; k++) {
//		int i, j;
//		double x;
//		in >> i >> j >> x;
//		i--; // C uses 0-indices, MatrixMarket uses 1-indices
//		j--;
//		if (x != 1) {
//			cerr << "Non-zero/one element in matrix." << endl;
//			exit(-1);
//		}
//		columns[j]->push_back(i);
//	}
//
//	// Sort all columns, just in case MatrixMarket file is corrupted
//	for (int j = 0; j < nCols; j++) {
//		std::sort(columns[j]->begin(), columns[j]->end());
//	}
//
//#ifdef DEBUG
//	cerr << "Read in sparse indicator matrix from " << fileName << endl;
//	cerr << "Spare matrix dimensions = " << nRows << " x " << nCols << endl;
//	cerr << "Number of non-zero elements = " << nEntries << endl;
//#endif
//
//}

CompressedDataMatrix::~CompressedDataMatrix() {
#ifdef DATA_AOS
	typedef std::vector<CompressedDataColumn*>::iterator CIterator;
	for (CIterator it = allColumns.begin(); it != allColumns.end(); ++it) {
		delete *it;
	}
#else
	typedef std::vector<real_vector*>::iterator RIterator;
	for (RIterator it = data.begin(); it != data.end(); ++it) {
		if (*it) {
			delete *it;
		}
	}

	typedef std::vector<int_vector*>::iterator IIterator;
	for (IIterator it = columns.begin(); it != columns.end(); ++it) {
		if (*it) {
			delete *it;
		}
	}
#endif
}

real CompressedDataMatrix::sumColumn(int column) {
#ifdef DATA_AOS
	cerr << "Not yet implemented.\n";
	exit(-1);
#else
	real sum = 0.0;
	if (getFormatType(column) == DENSE) {
		cerr << "Not yet implemented (DENSE)." << endl;
		exit(-1);
	} else if (getFormatType(column) == SPARSE) {
		cerr << "Not yet implemented (SPARSE)." << endl;
		exit(-1);
	} else { // is indiciator
		sum = columns[column]->size();
	}
	return sum;
#endif
}

void CompressedDataMatrix::printColumn(int column) {
#ifdef DATA_AOS
	cerr << "Not yet implemented.\n";
	exit(-1);
#else
	real_vector values;
	if (getFormatType(column) == DENSE) {
		values.assign(data[column]->begin(), data[column]->end());
	} else {
		bool isSparse = getFormatType(column) == SPARSE;
		values.assign(nRows, 0.0);
		int* indicators = getCompressedColumnVector(column);
		int n = getNumberOfEntries(column);
		for (int i = 0; i < n; ++i) {
			const int k = indicators[i];
			if (isSparse) {
				values[k] = data[column]->at(i);
			} else {
				values[k] = 1.0;
			}
		}
	}
	printVector(values.data(), values.size());
#endif
}

//template <class T>
//void CompressedDataMatrix::printVector(T values, const int size) {
//	cout << "[" << values[0];
//	for (int i = 1; i < size; ++i) {
//		cout << " " << values[i];
//	}
//	cout << "]" << endl;
//}



void CompressedDataMatrix::convertColumnToSparse(int column) {
#ifdef DATA_AOS
	cerr << "Not yet implemented.\n";
	exit(-1);
#else
	if (getFormatType(column) == SPARSE) {
		return;
	}
	if (getFormatType(column) == DENSE) {
		fprintf(stderr, "Format not yet support.\n");
		exit(-1);
	}

	while (data.size() <= column) {
		data.push_back(NULL);
	}
	if (data[column] == NULL) {
		data[column] = new real_vector();
	}

#if 1
	const real value = 1.0;
#else
	const real value = 2.0;
#endif

	data[column]->assign(nRows, value);
	formatType[column] = SPARSE;
#endif
}

void CompressedDataMatrix::convertColumnToDense(int column) {
#ifdef DATA_AOS
	cerr << "Not yet implemented.\n";
	exit(-1);
#else
	if (getFormatType(column) == DENSE) {
		return;
	}
	if (getFormatType(column) == SPARSE) {
		fprintf(stderr, "Format not yet support.\n");
		exit(-1);
	}

	while (data.size() <= column) {
		data.push_back(NULL);
	}
	if (data[column] == NULL) {
		data[column] = new real_vector();
	}
	data[column]->resize(nRows, static_cast<real>(0));

	int* indicators = getCompressedColumnVector(column);
	int n = getNumberOfEntries(column);
//	int nonzero = 0;
	for (int i = 0; i < n; ++i) {
		const int k = indicators[i];
//		cerr << " " << k;
//		nonzero++;

#if 1
		const real value = 1.0;
#else
		const real value = 2.0;
#endif
		data[column]->at(k) = value;
	}
//	cerr << endl;
//	cerr << "Non-zero count: " << nonzero << endl;
//	exit(0);
	formatType[column] = DENSE;
	delete columns[column]; columns[column] = NULL;
#endif
}

int CompressedDataMatrix::getNumberOfRows(void) const {
	return nRows;
}

int CompressedDataMatrix::getNumberOfColumns(void) const {
	return nCols;
}

void CompressedDataMatrix::setNumberOfColumns(int nColumns) {
	nCols = nColumns;
}

int CompressedDataMatrix::getNumberOfEntries(int column) const {
#ifdef DATA_AOS
	return allColcolumns[column]->getNumberOfEntries();
#else
	return columns[column]->size();
#endif
}

int* CompressedDataMatrix::getCompressedColumnVector(int column) const {
#ifdef DATA_AOS
	return allColumns[column]->getColumns();
#else
	return const_cast<int*>(&(columns[column]->at(0)));
#endif
}

void CompressedDataMatrix::removeFromColumnVector(int column, int_vector removeEntries) const{
	int lastit = 0;
	int_vector::iterator it1 = removeEntries.begin();
	int_vector::iterator it2 = columns[column]->begin();
	while(it1 < removeEntries.end() && it2 < columns[column]->end()){
		if(*it1 < *it2)
			it1++;
		else if(*it2 < *it1){
			it2++;
		}
		else{
			columns[column]->erase(it2);
			it2 = columns[column]->begin() + lastit;
		}
	}
}

void CompressedDataMatrix::addToColumnVector(int column, int_vector addEntries) const{
	int lastit = 0;
	int p = columns[column]->size() + addEntries.size();
	for(int i = 0; i < (int)addEntries.size(); i++)
	{
		int_vector::iterator it = columns[column]->begin() + lastit;
		if(columns[column]->size() > 0){
			while(*it < addEntries[i]){
				it++;
				lastit++;
			}
		}
		columns[column]->insert(it,addEntries[i]);
	}
}

real* CompressedDataMatrix::getDataVector(int column) const {
#ifdef DATA_AOS
	return allColumns[column]->getData();
#else
	return const_cast<real*>(data[column]->data());
#endif
}

void CompressedDataMatrix::getDataRow(int row, real* x) const {
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
}

//void CompressedDataMatrix::allocateMemory(int nCols) {
//	// Allocate some memory
////	columns = std::vector<int_vector*>(nCols);
//	columns.resize(nCols);
//	for (int j = 0; j < nCols; j++) {
////		columns[j] = int_vector(); // Create empty list
//		columns[j] = new int_vector();
//	}
//}

FormatType CompressedDataMatrix::getFormatType(int column) const {
#ifdef DATA_AOS
	allColumns[column]->getFormatType();
#else
	return formatType[column];
#endif
}

CompressedDataMatrix* CompressedDataMatrix::transpose(){
	CompressedDataMatrix* matTranspose = new CompressedDataMatrix();

	matTranspose->nRows = this->getNumberOfColumns();
	matTranspose->nCols = this->getNumberOfRows();

	bool flagDense = false;
	bool flagIndicator = false;
	for(int i = 0; i < nCols; i++){
		if(formatType[i] == DENSE)
			flagDense = true;
		if(formatType[i] == INDICATOR)
			flagIndicator = true;
	}
	matTranspose->formatType.resize(matTranspose->nCols);
	if(flagIndicator){
		for (int k = 0; k < matTranspose->nCols; k++) {
			int_vector* thisColumn = new int_vector();
			matTranspose->columns.push_back(thisColumn);
			matTranspose->formatType[k] = INDICATOR;
		}
	}

	if(flagDense){
		for (int k = 0; k < matTranspose->nCols; k++) {
			real_vector* thisData = new real_vector();
			matTranspose->data.push_back(thisData);
			matTranspose->formatType[k] = DENSE;
		}
	}

	if(flagIndicator && flagDense){
		for(int k = 0; k < matTranspose->nCols; k++){
			matTranspose->formatType[k] = SPARSE;
		}
	}


	for (int i = 0; i < matTranspose->nRows; i++) {
		if(formatType[i] == INDICATOR || formatType[i] == SPARSE){
			int rows = this->getNumberOfEntries(i);
			for (int j = 0; j < rows; j++) {
				matTranspose->columns[this->getCompressedColumnVector(i)[j]]->push_back(i);
				if(flagDense){
					if(formatType[i] == SPARSE)
						matTranspose->data[this->getCompressedColumnVector(i)[j]]->push_back(this->data[i]->at(j));
					else
						matTranspose->data[this->getCompressedColumnVector(i)[j]]->push_back(1.0);
				}
			}
		}
		else{
			for (int j = 0; j < nRows; j++) {
				if(flagIndicator){
					matTranspose->columns[j]->push_back(i);
				}
				matTranspose->data[j]->push_back(this->data[i]->at(j));
			}
		}
	}


	return matTranspose;
}