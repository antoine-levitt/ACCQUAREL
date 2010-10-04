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
#ifndef _TENSOR_NON_SPARSE_
#define _TENSOR_NON_SPARSE_

#include <iostream>
#include "tensor.h"

using namespace std;

class tensorNonSparse : public tensor
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
	 * The Method that computes the storage size for the tensor.
	 */
	virtual int computeStorageSize(void) const = 0;
	
	/**
	 * The Copy method.
	 */
	void copy(const tensorNonSparse & m);
	
	/**
	 * Method that frees the memory of the sorage.
	 */
	void freeStorage(void);
	
	/**
	 * The hash method.
	 *
	 * This method gives the position in the storage array of the
	 * element row x column of the tensor.
	 */
	virtual int hash(const int & row , const int & column , const int & width , const int & height) const = 0;

	/**
	 * Method SET to set a element of the storage array directly.
	 */
	void setStorage(const int & i , const double & value);
	
	/**
		* Method SET for the storage size.
	 */
	void setStorageSize(const int & storageSize);
	
public:
		
	/**
	 * The constructor.
	 */
	tensorNonSparse(void);
	
	/**
	 * The destructor.
	 */
	virtual ~tensorNonSparse(void);
	
	/**
	 * The Clear method.
	 */
	void clear(void);

	/**
	 * Method GET for the coefficient.
	 */
	virtual const double & getCoefficient(const int & row , const int & column , const int & width , const int & height) const;
	
	/**
	 * Method GET for the coefficient.
	 */
	double & getCoefficient(const int & row , const int & column , const int & width , const int & height);

	/**
	 * Method GET to access directly to the storage array.
	 */
	const double & getStorage(const int & i) const;
	
	/**
	 * Method GET to access directly to the storage array.
	 */
	double & getStorage(const int & i);
	
	/**
	 * Method GET for the Storage Size.
	 */
	const int & getStorageSize(void) const;

	/**
	 * Constantn operator to access the tensor as a storage array.
	 *
	 * @param item the position in the storage array of the element that should be acceded.
	 *
	 * @warning the item parameter shall be positive and strictly inferior to the storage size
	 * of the tensor non sparse.
	 */
	const double & operator[](const int & item) const;

	/**
	 * Operator to access the tensor as a storage array.
	 *
	 * @param item the position in the storage array of the element that should be acceded.
	 *
	 * @warning the item parameter shall be positive and strictly inferior to the storage size
	 * of the tensor non sparse.
	 */
	double & operator[](const int & item);

	/**
	 * Constant operator to access the tensor as a storage array.
	 *
	 * @param row the row.
	 * @param column the column.
	 * @param width the width.
	 * @param height the height.
	 *
	 * @warning all the parameters shall be in the defintion range.
	 */
	const double & operator()(const int & row , const int & column , const int & width , const int & height) const;

	/**
	 * Operator to access the tensor as a storage array.
	 *
	 * @param row the row.
	 * @param column the column.
	 * @param width the width.
	 * @param height the height.
	 *
	 * @warning all the parameters shall be in the defintion range.
	 */
	double & operator()(const int & row , const int & column , const int & width , const int & height);

	/**
	 * Method SET : this methos sets all the coefficient of the 
	 * tensor to a particular value.
	 */
	void setAllCoefficients(const double & value);
	
	/**
	 * Method SET for the coefficient.
	 */
	void setCoefficient(const int & row , const int & column , const int & width , const int & height , const double & value);
	
	/**
	 * Method SET for the Size of the tensor.
	 */
	void setTensorSize(const int & nbrOfRows , const int & nbrOfColumns , const int & nbrOfWidth , const int & nbrOfHeights);
};

#endif
