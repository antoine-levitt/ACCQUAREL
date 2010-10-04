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
#ifndef _VECTOR_
#define _VECTOR_

#include <iostream>
#include<math.h>
#include "matrix.h"
#include "matrixFull.h"
#include "matrixSymetric.h"
using namespace std;

class vector
{

private:
	
	/**
	 * The storage array for data.
	 */
	double * Coefficients;
	
	/**
	 * The size of the vector.
	 */
	int NbrOfRows;
	
protected:
	
	/**
	 * this method performs memory allocation to store 
	 * the good nbr of rows.
	 */
	void allocStorage(const int & nbrOfRows);
	
	/**
	 * This method perform the copy of an object vector.
	 */
	void copy(const vector & v);
	
	/**
	 * This method frees the memory of the vector.
	 */
	void freeStorage(void);
	
public:
	
	/**
	 * Default constructor.
	 */
	vector(void);
	
	/**
	 * Constructor with size specification.
	 */
	vector(int nbrOfRows);

	/**
	 * Copy constructor.
	 */
	vector(const vector & v);

	/**
	 * destructor.
	 */
	virtual ~vector(void);

	/**
	 * This method clears the vector.
	 */
	void clear(void);
	
	/**
	 * This method returns true if the vector is empty
	 * false if not.
	 */
	bool empty(void) const;

	/**
	 * Method GET for the coefficient of the vector.
	 */
	const double & getCoefficient(const int & row) const;
	
	/**
	 * Method GET for the coefficient of the vector.
	 */
	double & getCoefficient(const int & row);

	/**
	 * Method GET for the umber of rows.
	 */
	int getNbrOfRows(void) const;
	
	/**
	 * Operator to acces coefficients of the vector.
	 */
	const double & operator[] (const int & row) const;

	/**
	 * Operator to acces coefficients of the vector.
	 */
	double & operator[] (const int & row);

	/**
	 * Method SET for the coefficient of the vector.
	 */
	void setCoefficient(const int & row , const double & value);
	
	/**
	 * Method SET for the number of rows.
	 */
	void setNbrOfRows(const int & nbrOfRows);

	/**
	 * Method SET for the vector size.
	 *
	 * This is exactly the same thing as calling setNbrOfRows(void).
	 */
	void setVectorSize(int size);
	
	/**
	 * Method that computes the norm
	 */
    double Norm(void) const; 
	
	/**
    * The friend method that computes the scalar product of two vectors
	*/
	friend const double vectorScalarProduct(const vector & u,const vector & v); 
	
	/**
    * The friend method that computes the L2 scalar product of two vectors with respect to an overlap matrix
	*/
	 friend const double L2vectorScalarProduct(const vector & u,const vector & v, const matrixSymetric & Moverlap); 
	 
	 /**
    * The friend method that computes the L2 norm with respect to an overlap matrix
	*/
	 friend const double L2vectorNorm(const vector & u, const matrixSymetric & Moverlap); 

	/**
	 * The friend method that computes the product of a matrix and a vector
	 */
	 friend const vector ProductMatrixVector(const matrixFull & A, const vector & v);
	 
	 /**
	 * The friend method that solves the linear equation Ax=LUx=b
	 */
	 friend const void LUsolver(vector & x,const vector & b, const matrixFull & A);
	 
	 /**
	 * The friend method that computes the inverse of a matrix from its LU decomposition
	 */
	 friend const void LUinverse( matrixFull & inv, const matrixFull & A);
	 
	 /**
	 * The method that redefine the operator =
	 */
	 vector operator = (const vector & v);
	 
	/**
	 * Method that write a vector in a stream
	 */
	void writeHuman(ostream & outStream , int width = 8 , int precision = 4) const;
};


inline ostream & operator<< (ostream & outStream , const vector & v)
{
	v.writeHuman(outStream);
	return outStream;
}
#endif

