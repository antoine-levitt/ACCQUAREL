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
#ifndef _MATRIX_SYMETRIC_
#define _MATRIX_SYMETRIC_

#include <assert.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "matrixNonSparse.h"
#include "tensorSymetric.h"

using namespace std;

class matrixSymetric : public matrixNonSparse
{
private:

protected:

/**
 * Method that computes the storage size.
 */
virtual int computeStorageSize(void) const;

/**
* Method that computes the storage size.
 */
virtual int hash(const int & row , const int & column ) const;

public:

/**
 * The default constructor.
 */
matrixSymetric(void);

/**
 * The constructor with size specification.
 *
 * This method constructs a symetric matrix that will be aible to
 * store nbrOfRows x nbrOfRows elements.
 */
matrixSymetric(int nbrOfRows);

/**
 * The constructor with size specification.
 *
 * This method constructs a symetric matrix that will be aible to
 * store nbrOfRows x nbrOfRows elements.
 */
matrixSymetric(const matrixSymetric & matrixSym);

/**
 * The destructor.
 */
virtual ~matrixSymetric(void);	

	/**
	 * Affectation operator.
	 *
	 * @param matrixSym the matrix symetric to add to the current
	 * object.
	 */
	matrixSymetric & operator= (const matrixSymetric & matrixSym);

	/**
	 * Operator to acces the unary addition of two matrixes.
	 *
	 * @param matrixSym the matrix symetric to add to the current
	 * object.
	 *
	 * @warning sizes of the calling object and matrixSym must match.
	 */
	matrixSymetric & operator+= (const matrixSymetric & matrixSym);

	/**
	 * Operator to acces the binary addition of two matrixes.
	 *
	 * @param matrixSym the matrix symetric to add to the current
	 * object.
	 *
	 * @warning sizes of the calling object and matrixSym must match.
	 */
	matrixSymetric operator+ (const matrixSymetric & matrixSym) const;

	/**
	 * Operator to acces the unary soustraction of two matrixes.
	 *
	 * @param matrixSym the matrix symetric to soustract to the current
	 * object.
	 *
	 * @warning sizes of the calling object and matrixSym must match.
	 */
	matrixSymetric & operator-= (const matrixSymetric & matrixSym);

	/**
	 * Operator to acces the binary soustraction of two matrixes.
	 *
	 * @param matrixSym the matrix symetric to soustract to the current
	 * object.
	 *
	 * @warning sizes of the calling object and matrixSym must match.
	 */
	matrixSymetric operator- (const matrixSymetric & matrixSym) const;

	/**
	 * Operator to acces the unary multiplication of a matrix and a scalar.
	 *
	 * @param scalar the scalar value the matrix shall be multiplied with.
	 */
	matrixSymetric & operator*= (const double & scalar);

	/**
	 * Operator to acces the binary multiplication of a matrix and a scalar.
	 *
	 * @param scalar the scalar value the matrix shall be multiplied with.
	 */
	matrixSymetric operator* (const double & scalar) const;

		/**
	 * Operator to acces the unary division of a matrix and a scalar.
	 *
	 * @param scalar the scalar value the matrix shall be divided with.
	 */
	matrixSymetric & operator/= (const double & scalar);

	/**
	 * Operator to acces the binary division of a matrix and a scalar.
	 *
	 * @param scalar the scalar value the matrix shall be divided with.
	 */
	matrixSymetric operator/ (const double & scalar) const;

/**
 * Method read : This method reads a matrix in
 * an unformated file.
 */
void readText(const string & fileName);

 /**
  * Method read : This method reads a matrix in
  * an unformated file.
  */
 void readText(istream & inStream);

 /**
	* Method SET for the matrix Size.
	* 
	*
	* @param the number of rows (equal to the number of columns) of
	* that the object will be able to contain.
	*/
void setMatrixSize(const int & nbrOfRows);

 /**
	* Method SET for the matrix Size.
	*
	* This method is here to reimplement the method from the 
	* matrixNonSparse classe. If the number of rows and the number
	* of columns differs this method will fail.
	* 
	* @param nbrOfRows the number of rows of the matrix.
	*
	* @param nbrOfColumns the number of columns of the matrix.
	*
	* @warning this will fail if the numbers of rows is different from the number of columns.
	*/
 void setMatrixSize(const int & nbrOfRows , const int & nbrOfColumns);

 /**
	* Method that writes a matrix symetric in a stream.
	*/
 void writePretty(ostream & outStream , int width = 10 , int precision = 6) const;

 /**
	* Method that writes the symetric Matrix in a file for that can be 
	* read with the method readText.
	*/
 void writeText(ostream & outStream , int precision = 6) const;
};

/**
 * extern operator to multiply a matrix symetric and a scalar
 */
extern matrixSymetric operator * (const double & scalar , const matrixSymetric & matrixSym);

#endif

