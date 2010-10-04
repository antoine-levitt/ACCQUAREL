/* 
* The martix library of the A.S.P.I.C. 
 * Written and directed by François Lodier francois.lodier@gmail.com.
 * Modified by Frederic Legoll nov 06
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
#ifndef _MATRIX_FULL_
#define _MATRIX_FULL_

#include <iomanip>
#include <iostream>
#include <math.h>
#include "matrixNonSparse.h"

using namespace std;

class matrixFull : public matrixNonSparse
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
	virtual int hash( const int & row , const int & column ) const;

public:
	
	/**
	 * The constructor.
	 */
	matrixFull(void);

	/**
	 * The constructor.
	 */
	matrixFull(const int & nbrOfRows , const int & nbrOfColumns);
	
		/**
	 * The constructor.
	 */
	matrixFull(const matrixFull & mFull);
	
	/**
	 * The destructor.
	 */
	virtual ~matrixFull(void);	
	
	/**
		* Méthode bête pour faire la multiplication de deux martrices.
	 *
	 * Cette méthode permet de Multiplier la matrice A avec la 
	 * Matrice B , à partir du moment où les tailles correspindent.
	 * Le résultat est stocké dans la matrices qui appelle la méthode
	 */
	void multiply(const matrix & A , const matrix & B);
	
	/**
	 * The method that computes the norm of a column of a matrix
	 */
     friend double MatrixFullNormColumn(const matrixFull & A, int column);
	 
	 /**
	 * The method that computes the scalar product of two columns of a matrix
	 */
     friend double MatrixFullScalarProduct(const matrixFull & A, int i, int j);
	
	/**
	 * The method to write matrixes in the prettiest way.
	 */
	void writePretty(ostream & outStream , int width = 10 , int precision = 6) const;
	
	/**
	 * The method that takes the transposition of a full matrix
	 */
	 matrixFull TransposeMatrix(void) const;
	 
	 /**
	 * The method that takes the signature of an element of a matrix 
	 */
	 int MatrixSignElement(int i, int j) const;
	 
	 /**
	 * The method that finds a minor of a square matrix 
	 */
	 void MinorForSquareMatrix(matrixFull & M1, int line, int column) const;
	 
	 /**
	 * The method that computes the determinant of a matrix 
	 */
	 double Determinant(void) const;
	 
	 /**
	 * The method that computes the comatrix of a matrix 
	 */
	 matrixFull Comatrix(void) const;
	 
	 /**
	 * The method that computes the inverse of a full matrix
	 */
	 matrixFull Inverse(void) const;
	 
	 /**
	 * The method that computes the LU decomposition of a matrix
	 */
	 void LUdecomposition(matrixFull & LU, int *indx, double d) const;

	/**
	 * Affectation operator.
	 */
	matrixFull & operator= (const matrixFull & mFull);
};


























#endif

