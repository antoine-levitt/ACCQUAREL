/* 
 * The martix library of the A.S.P.I.C. 
 * Written and directed by François Lodier francois.lodier@gmail.com.
 *
 * Copyright (C) 2005  François Lodier
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include "generalizedEigenSolver.h"

//////////////////////////////////////////////////////////////////////////////////////////////
// This is the declaration of the external lapack routines that solves the genralized problem.
//////////////////////////////////////////////////////////////////////////////////////////////
extern "C" {
int dspgv_(
					 int *itype, 
					 char *jobz, 
					 char *uplo, 
					 int * n, 
					 double *ap, 
					 double *bp, 
					 double *w, 
					 double *z__, 
					 int *ldz, 
					 double *work, 
					 int *info);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// The memory allocation method.
/////////////////////////////////////////////////////////////////////////////////////////////////
generalizedEigenSolver::~generalizedEigenSolver(void) 
{
	if(!empty()) {
		free();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// The memory allocation method.
/////////////////////////////////////////////////////////////////////////////////////////////////
void generalizedEigenSolver::alloc(const int & nbrOfRows) 
{
	assert(nbrOfRows > 0);
	assert(empty());

	EigenVectors.setSizes(nbrOfRows);
	
	for(int i = 0 ; i < nbrOfRows ; i++)
		EigenVectors[i].setVectorSize(nbrOfRows);

	EigenValues.setSizes(nbrOfRows);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// The clear method.
/////////////////////////////////////////////////////////////////////////////////////////////////
void generalizedEigenSolver::clear(void) 
{
	if(empty()) {
		return;
	} else {
		free();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// The empty method.
/////////////////////////////////////////////////////////////////////////////////////////////////
bool generalizedEigenSolver::empty(void) const
{
	if(getNbrOfEigenValues() == 0) {
		return true;
	} else {
		return false;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// The memory allocation method.
/////////////////////////////////////////////////////////////////////////////////////////////////
void generalizedEigenSolver::free(void) 
{
	if(empty()) {
		return;
	}

	EigenVectors.clear();
	EigenValues.clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// The acces to the eigenvalues.
/////////////////////////////////////////////////////////////////////////////////////////////////
const double & generalizedEigenSolver::getEigenValue(const int & item) const
{
	assert(item < getNbrOfEigenValues());
	
	return EigenValues[item];
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// The acces to the eigenvalues.
/////////////////////////////////////////////////////////////////////////////////////////////////
const containor<double> & generalizedEigenSolver::getEigenValues(void) const
{
	return EigenValues;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// The acces to the eigenvector.
/////////////////////////////////////////////////////////////////////////////////////////////////
const vector & generalizedEigenSolver::getEigenVector(const int & item) const
{
	assert(item < getNbrOfEigenValues());
	
	return EigenVectors[item];
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// The acces to the size of the problem.
/////////////////////////////////////////////////////////////////////////////////////////////////
const containor<vector> & generalizedEigenSolver::getEigenVectors(void) const
{
	return EigenVectors;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// The acces to the size of the problem.
/////////////////////////////////////////////////////////////////////////////////////////////////
const int & generalizedEigenSolver::getNbrOfEigenValues(void) const
{
	return EigenVectors.getNbrOfRows();
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Method that solves the generalized eigen values problem Ax = lambda.Bx
//////////////////////////////////////////////////////////////////////////////////////////////
void generalizedEigenSolver::solve(const matrixSymetric & A , const matrixSymetric & B) 
{
	assert(A.getNbrOfRows() == B.getNbrOfRows());


	bool storeResults;
	int job , info , nbrOfRows = A.getNbrOfRows() , i , j;
	double * arrayA , * arrayB , * eigenValues , *eigenVectors , *work;
	char v , u;
	
	if(!empty()) {
		clear();
	}

	if(nbrOfRows == 0) {
		return;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialisation of the lapack datas.
	arrayA = new double [A.getStorageSize()];
	for(i=0 ; i < A.getStorageSize() ; i++) 
		arrayA[i] = A[i];
	
	arrayB = new double [B.getStorageSize()];
	for(i=0 ; i < B.getStorageSize() ; i++) 
		arrayB[i] = B[i];

	eigenValues = new double[nbrOfRows];
	eigenVectors = new double[nbrOfRows * nbrOfRows];
	work = new double [3* nbrOfRows];
	
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Calling Lapack.
	job=1;
	v='V';
	u='U';
	dspgv_(&job,&v,&u,&nbrOfRows,arrayA, arrayB , eigenValues , eigenVectors ,&nbrOfRows , work , &info);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Catching errors.
	switch(info) {
		case 0 :
			storeResults = true;
			break;

		default:

			cerr << "Error : in void generalizedEigenSolver::solve(const matrixSymetric & A , const matrixSymetric & B)" << endl;
			cerr << "Error : Unknown Error" << endl;
			storeResults = false;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Storing results if needed.
	if(storeResults) {
		alloc(nbrOfRows);

		for(i=0 ; i < nbrOfRows ; i++) {
			EigenValues[i] = eigenValues[i];
		}
	
		for(i=0 ; i < nbrOfRows ; i++) {
			for(j=0 ; j < nbrOfRows ; j++) {
				EigenVectors[j][i] = eigenVectors[j * nbrOfRows + i];
			}
		}

	}

	////////////////////////////////////////////////////////////////////////////////////
	// Here we delete all the arrays.
	////////////////////////////////////////////////////////////////////////////////////
	delete [] arrayA;
	delete [] arrayB;
	delete [] eigenValues;
	delete [] eigenVectors;
	delete [] work;
}
