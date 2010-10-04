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
#ifndef _MATRIX_NON_SPARSE_
#define _MATRIX_NON_SPARSE_

#include <iostream>
#include "matrix.h"

using namespace std;

class matrixNonSparse : public matrix
{
	
private:
	
	/**
	 * The array containing the coefficients.
	 */
	double * Storage;
	
	/**
	 * The Size of the storage.
	 */
	int StorageSize;

protected:

	/**
	 * Method that allocates memory for the given storage size.
	 */
	void allocStorage(const int & storagesize);
		
	/**
	 * The Method that computes the storage size for the matrix.
	 */
	virtual int computeStorageSize(void) const = 0;
	
	/**
	 * The Copy method.
	 */
	void copy(const matrixNonSparse & m);
	
	/**
	 * Method that frees the memory of the sorage.
	 */
	void freeStorage(void);
	
	/**
	 * The hash method.
	 *
	 * This method gives the position in the storage array of the
	 * element row x column of the matrix.
	 */
	virtual int hash(const int & row , const int & column) const = 0;
	
public:
		
	/**
	 * The constructor.
	 */
	matrixNonSparse(void);
	
	/**
	 * The destructor.
	 */
	virtual ~matrixNonSparse(void);
	
	/**
	 * The method to add two matrixes.
	 *
 * This method adds matrixA and matrixB and stores the results
 * in the calling object. The correct memory allocation is performed
 * on the current object to store the result.
 *
 * @warning sizes of the matrixA and matrixB must match.
 *
 * @warning matrix types of matrixA and matrixB must match.
 *
 * @warning this method is a ! DUMMY ! one : User must ensure
 * that matrixA, matrixB and the calling aobject have the same matrix
 * type. If types mismatch an assert will abort the program 
 * in debug mode and results may be surprising in the release mode.
 */
	virtual void add(const matrixNonSparse & matrixA , const matrixNonSparse & matrixB);

/**
 *
 * The to unary add two symetric matrixes.
 *
 * This method adds matrixA and the calling object. The result is stored
 * in the calling object. User must be aware, that this is the ! DUMMY !
 * parent method.
 *
 * @warning sizes of the matrixA and object must match to perform the sum. 
 *
 * @warning this method is a ! DUMMY ! one : User must ensure
 * that matrixA and the calling aobject have the same matrix
 * type. If types mismatch an assert will abort the program 
 * in debug mode and results may be surprising in the release mode.
 */
void add(const matrixNonSparse & matrixA);

	/**
	 * The Clear method.
	 */
	void clear(void);

	/**
   * Binary matrix - scalar division.
   * This method divides the matrixA with the scalar. The result
	 * of the operation is stored in the calling object.
	 *
	 * @param matrixA the matrix to be divided by the scalar value.
	 * 
	 * @param scalar the scalar with which the matrixA shall be divided.
	 *
	 * @warning the scalar value shall not be null.
	 *
	 * @warning this method is a dummy one : User must ensure
	 * that matrixA and the calling aobject have the same matrix
	 * type. If types mismatch an assert will abort the program 
	 * in debug mode and results may be surprising in the release mode.
   */
	virtual void divide(const matrixNonSparse & matrixA , const double & scalar);
	
	/**
   * Unary scalar - matrix division.
   * This method divides the calling object with the scalar. The result
	 * of the operation is stored in the calling object.
   *
	 * @param scalar the scalar to divide the calling object with.
	 *
	 * @warning the scalar value shall not be null.
	 */
	virtual void divide(const double & scalar);

	/**
	 * Method GET for the coefficient.
	 */
	virtual const double & getCoefficient(const int & row , const int & column) const;

	/**
	 * Method GET for the coefficient.
	 */
	virtual double & getCoefficient(const int & row , const int & column);
	
	/**
	 * Method GET for the coefficient.vv
	 */
	const double * getStorage(void) const;
	
	/**
	 * Method GET to access directly to the storage array.
	 */
	const double & getStorage(const int & i) const;
	
	/**
	 * Method GET to access directly to the storage array.
	 *
	 * @param item the position in the storage array of the element to find.
	 */
	double & getStorage(const int & item);
	
	/**
	 * Method GET for the Storage Size.
	 *
	 * @return the size allocated for the storage.
	 */
	const int & getStorageSize(void) const;
	
	/**
   * The method to a matrix and a scalar.
   * This method multiply the matrixA with the scalar. The result
	 * of the operation is stored in the calling object.
	 *
	 * @warning this method is a dummy one : User must ensure
	 * that matrixA and the calling aobject have the same matrix
	 * type. If types mismatch an assert will abort the program 
	 * in debug mode and results may be surprising in the release mode.
   */
	virtual void multiply(const matrixNonSparse & matrixA , const double & scalar);

	/**
   * The method to a matrix and a scalar.
   *
   * This method multiply the calling object with the scalar.
   */
	virtual void multiply(const double & scalar);

	/**
	 * Constant access operator for the elements of the matrix as a
	 * array.
	 *
	 * @param item the position in the storage array of the element that 
	 * shall be acceded.
	 *
	 * @return a constant reference for the item element stored.
	 *
	 * @warning all parameters shall be in the appropriate range.
	 */
	const double & operator[] (const int & item) const;
	
	/**
	 * Access operator for the elements of the matrix as a
	 * array.
	 *
	 * @param item the position in the storage array of the element that 
	 * shall be acceded.
	 *
	 * @return a constant reference for the item element stored.
	 *
	 * @warning all parameters shall be in the appropriate range.
	 */
	double & operator[] (const int & item);
	
	/**
	 * Constant access operator for the elements of the matrix.
	 *
	 * @param row the row of the element to access.
	 * @param column the column of the element to access.
	 *
	 * @return a constant reference for the element row x column of the matrix.
	 *
	 * @warning all parameters shall be in the appropriate range.
	 */
	const double & operator() (const int & row , const int & column) const;
	
	/**
	 * Access operator for the elements of the matrix.
	 *
	 * @param row the row of the element to access.
	 *
	 * @param column the column of the element to access.
	 */
	double & operator() (const int & row , const int & column);
		
	/**
	 * Method SET for the coefficient.
	 */
	void setCoefficient(const int & row , const int & column , const double & value);
	
  /**
	 * Method SET : this methos sets all the coefficient of the 
	 * matrix to a particular value.
	 *
	 * @param value the value that all the coefficients from the 
	 * matrix will have.
	 */
	void setCoefficients(const double & value);
	
	/**
	 * Method SET for the Size of the matrix.
	 *
	 * @param nbrOfRows the number of rows that the matrix will be abble to contain
	 * @param nbrOfColumns the number of columns that the matrix will be abble to contain
	 */
	void setMatrixSize(const int & nbrOfRows , const int & nbrOfColumns);
	
	/**
	 * Method SET to set a element of the storage array directly.
	 *
	 * @param item the position in the storage array of the element
	 * the value will be changed.
	 *
	 * @param value the value of that the element will take.
	 */
	void setStorage(const int & i , const double & value);
	
	/**
	 * Method SET for the storage size.
	 */
	void setStorageSize(const int & storageSize);

/**
 * Binary soustraction of two matrixes.
 *
 * This method soustracts matrixA and matrixB and stores the results
 * in the calling object. The correct memory allocation is performed
 * on the current object to store the result.
 *
 * @warning sizes of the matrixA and matrixB must match.
 *
 * @warning matrix types of matrixA and matrixB must match.
 *
 * @warning this method is a ! DUMMY ! one : User must ensure
 * that matrixA, matrixB and the calling aobject have the same matrix
 * type. If types mismatch an assert will abort the program 
 * in debug mode and results may be surprising in the release mode.
 */
	virtual void soustract(const matrixNonSparse & matrixA , const matrixNonSparse & matrixB);

/**
 * Unary soustraction of matrixes.
 *
 * This method soustracts matrixA to the calling object. The result is stored
 * in the calling object. User must be aware, that this is the ! DUMMY !
 * parent method.
 *
 * @warning sizes of the matrixA and object must match to 
 * perform the soustraction. 
 *
 * @warning this method is a ! DUMMY ! one : User must ensure
 * that matrixA and the calling aobject have the same matrix
 * type. If types mismatch an assert will abort the program 
 * in debug mode and results may be surprising in the release mode.
 */
void soustract(const matrixNonSparse & matrixA);


};

#endif

