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
#ifndef _MATRIX_
#define _MATRIX_

#include <assert.h>
#include <iostream>
#include <string>
using namespace std;

/**
 * The front class for all the matrix.
 */
class matrix
{

public:
	
	
	/**
	 * The types of matrixes.
	 */
	typedef enum _matrixType {
		FULL,
		NONE,
		SPARSE,
		SYMETRIC_FULL,
		SYMETRIC_SPARSE
} matrixType;


/**
* Function that converts a matrix type to a string.
 */
static const string matrixType2String(const matrix::matrixType & type);

/**
* Function that converts a string to a matrix type.
 */
static const matrix::matrixType string2MatrixType(const string & matrixTypeStr);


private:
		
	/**
	 * The number of columns in the matrix.
	 */
	int NbrOfColumns;
	
	/**
	 * The number of rows in the matrix.
	 */
	int NbrOfRows;
	
	/**
	 * The Type of the matrix.
	 */
	matrixType Type;
	
protected:

	/**
	 * The copy Method.
	 */
 void copy(const matrix & m);

	/**
	 * Method SET for the number of columns.
	 */
	void setNbrOfColumns(const int & nbrOfColumns);
	
	/**
	 * Method SET for the number of rows.
	 */
	void setNbrOfRows(const int & nbrOfRows);
	
	/**
	 * Method SET for the type of the matrix.
	 */
	void setMatrixType(const matrixType & type);
	
	
public :	
	
	
	/**
	 * The default constructor.
	 */
	matrix(void);
	
	/**
	 * The destructor.
	 */
	virtual ~matrix(void);

	/**
	 * Method Clear : this method makes the current object an empty one.
	 */
	void clear(void);
	
	/**
	 * Method to know if the matrix object is empty.
	 */
	bool empty(void) const;
	
	/**
	 * The GET methd for the coefficients of the matrix.
	 *
	 * This is a pure virtual method.
	 */
	virtual const double & getCoefficient(const int & row , const int & column)  const=0;
	
	/**
	 * Method GET for the number of rows.
	 */
	const int & getNbrOfRows(void) const;
	
	/**
	 * Method GET for the number of columns.
	 */
	const int & getNbrOfColumns(void) const;
	
	/**
	 * Method GET for the type of matrix.
	 */
	const matrixType & getMatrixType(void) const;

	/**
	 * The dummy method that computes the trace of the matrix.
	 *
	 * This method is virtual so it can be reimplemented in the 
	 * children class to increase efficiency.
	 *
	 * @return the trace of the matrix.
	 */
	virtual double trace(void) const;

};


/**
 * The operator to pretty print the matrix types.
 */
extern ostream & operator<<(ostream & outStream , const matrix::matrixType & type);

#endif

