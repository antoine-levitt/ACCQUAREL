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
#include "matrixSymetric.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The default constructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixSymetric::matrixSymetric(void)
: matrixNonSparse()
{
	setMatrixType(SYMETRIC_FULL);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The constructor with the sizes.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixSymetric::matrixSymetric(int nbrOfRows)
: matrixNonSparse()
{
	setMatrixType(SYMETRIC_FULL);
	setMatrixSize(nbrOfRows);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The constructor with the sizes.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixSymetric::matrixSymetric(const matrixSymetric & matrixSym)
: matrixNonSparse()
{
	setMatrixType(SYMETRIC_FULL);
	copy(matrixSym);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The destructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixSymetric::~matrixSymetric(void)
{
	;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Method that computes the total storage size.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int matrixSymetric::computeStorageSize(void) const
{
	assert(getNbrOfRows() == getNbrOfColumns());
	return getNbrOfRows() * (getNbrOfColumns() + 1) / 2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Method that computes the position in the storage array of the element row x column.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int matrixSymetric::hash(const int & row , const int & column) const
{
	assert(getNbrOfColumns() == getNbrOfColumns());
	
	assert( row >= 0 );
	assert( row < getNbrOfRows() );
	
	assert( column >= 0 );
	assert( column < getNbrOfColumns() );
	
	int orderRow , orderColumn;
	
	if( row > column ) {
		orderRow = column;
		orderColumn = row;
	} else {
		orderRow = row;
		orderColumn = column;
	}
	
	return  orderColumn * (orderColumn + 1) / 2 + orderRow;
}


/**
 * Affectation operator.
 *
 * @param matrixSym the matrix symetric to add to the current
 * object.
 */
matrixSymetric & matrixSymetric::operator= (const matrixSymetric & matrixSym)
{
	copy(matrixSym);
	return *this;
}

/**
* Operator to acces the unary addition of two matrixes.
*
* @param matrixSym the matrix symetric to add to the current
* object.
*
* @warning sizes of the calling object and matrixSym must match.
*/
matrixSymetric & matrixSymetric::operator += (const matrixSymetric & matrixSym)
{
	assert(matrixSym.getNbrOfRows() == getNbrOfRows());
	add(matrixSym);
	return *this;
}

/**
* Operator to acces the binary addition of two matrixes.
*
* @param matrixSym the matrix symetric to add to the current
* object.
*
* @warning sizes of the calling object and matrixSym must match.
*/
matrixSymetric matrixSymetric::operator+ (const matrixSymetric & matrixSym) const
{
	assert(matrixSym.getNbrOfRows() == getNbrOfRows());
	matrixSymetric tmp;
	tmp.add(*this,matrixSym);
	return tmp;
}

/**
* Operator to acces the unary soustraction of two matrixes.
*
* @param matrixSym the matrix symetric to soustract to the current
* object.
*
* @warning sizes of the calling object and matrixSym must match.
*/
matrixSymetric & matrixSymetric::operator-= (const matrixSymetric & matrixSym)
{
	assert(matrixSym.getNbrOfRows() == getNbrOfRows());
	soustract(matrixSym);
	return *this;
}

/**
* Operator to acces the binary soustraction of two matrixes.
*
* @param matrixSym the matrix symetric to soustract to the current
* object.
*
* @warning sizes of the calling object and matrixSym must match.
*/
matrixSymetric matrixSymetric::operator- (const matrixSymetric & matrixSym) const
{
	assert(matrixSym.getNbrOfRows() == getNbrOfRows());
	matrixSymetric tmp;
	tmp.soustract(*this,matrixSym);
	return tmp;
}

/**
* Operator to acces the unary multiplication of a matrix and a scalar.
*
* @param scalar the scalar value the matrix shall be multiplied with.
*/
matrixSymetric & matrixSymetric::operator*= (const double & scalar)
{
	multiply(scalar);
	return *this;
}

/**
* Operator to acces the binary multiplication of a matrix and a scalar.
*
* @param scalar the scalar value the matrix shall be multiplied with.
*/
matrixSymetric matrixSymetric::operator* (const double & scalar) const
{
	matrixSymetric tmp;
	tmp.multiply(*this,scalar);
	return tmp;
}

/**
 * extern operator to multiply a matrix symetric and a scalar
 */
matrixSymetric operator * (const double & scalar , const matrixSymetric & matrixSym)
{
	return matrixSym*scalar;
}

/**
* Operator to acces the unary division of a matrix and a scalar.
*
* @param scalar the scalar value the matrix shall be divided with.
*/
matrixSymetric & matrixSymetric::operator/= (const double & scalar)
{
	assert(scalar!= 0);
	divide(scalar);
	return *this;
}

/**
* Operator to acces the binary division of a matrix and a scalar.
*
* @param scalar the scalar value the matrix shall be divided with.
*/
matrixSymetric matrixSymetric::operator/ (const double & scalar) const
{
	assert(scalar!= 0);
	matrixSymetric tmp;
	tmp.divide(*this,scalar);
	return tmp;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method that reads a matrix in a file.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void matrixSymetric::readText(const string & fileName)
{
	ifstream fileStream;
	
	fileStream.open(fileName.c_str());
	
	if(!fileStream.is_open()) {
		cerr << "Error : in void matrixSymetric::readText(const string & fileName)." << endl;
		cerr << "Error : unable to open file with name " << fileName << endl;
		cerr << "Error : aborting reading now." << endl;
		clear();
		return;
	}
	
	readText(fileStream);
	
	if(fileStream.fail()) {
		cerr << "Error : in void matrixSymetric::readText(const string & fileName)." << endl;
		cerr << "Error : reading fails for file with name " << fileName << endl;
		clear();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that reads a matrix symetric in a stream.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void matrixSymetric::readText(istream & inStream) 
{

	int nbrOfRows;
	double value;
	
	inStream >> nbrOfRows;
	
	if(inStream.fail()) {
		cerr << "Error : in void matrixSymetric::readTxt(istream & inStream)." << endl;
		cerr << "Error : unable to read the number of rows in stream." << endl;
		cerr << "Error : aborting reading" << endl;
		clear();
		return;
	}
	
	setMatrixSize(nbrOfRows);
	
	for(int i=0 ; i<getStorageSize() ; i++) {
		inStream >> value;

		if(inStream.fail()) {
			cerr << "Error : in void matrixSymetric::readTxt(istream & inStream)." << endl;
			cerr << "Error : unable to read data number " << i << "." << endl;
			cerr << "Error : aborting reading" << endl;
			clear();
			return;
		}
		
		setStorage(i,value);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that sets the size of the matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void matrixSymetric::setMatrixSize(const int & nbrOfRows) 
{
	matrixNonSparse::setMatrixSize(nbrOfRows , nbrOfRows);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that sets the size of the matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void matrixSymetric::setMatrixSize(const int & nbrOfRows , const int & nbrOfColumns)
{
	assert( nbrOfRows == nbrOfColumns );
	matrixNonSparse::setMatrixSize(nbrOfRows , nbrOfColumns);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that writes down the matrix in a human readable way
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void matrixSymetric::writePretty(ostream & outStream , int width , int precision) const
{
	int i , j , nbrOfRows;
	

	nbrOfRows = getNbrOfRows();
	
	for(i=0 ; i < nbrOfRows ; i++) {
		
		for(j=0 ; j < nbrOfRows ; j++) {
		
			if(j < i)
				outStream << setw(width) << "x";
			else 
				outStream << setw(width) << setprecision(precision) << getCoefficient(i,j);
		
		}
		
		outStream << endl;
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that writes the matrix with the only data that it contains.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void matrixSymetric::writeText(ostream & outStream , int precision) const
{
	outStream << getNbrOfRows() << endl << endl;

	for(int i =0 ; i < getStorageSize() ; i++) {
		outStream << setprecision(precision) << getStorage(i) << endl;
	}
}

