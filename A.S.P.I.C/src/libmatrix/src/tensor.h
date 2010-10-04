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
#ifndef _TENSOR_
#define _TENSOR_

#include <assert.h>
#include <iostream>
#include <string>
using namespace std;


/**
 * The front class for all the tensor.
 */
class tensor
{

public :
	typedef enum _tensorType
{
	FULL,
	NONE,
	SPARSE,
	SYMETRIC_FULL,
	SYMETRIC_SPARSE
} tensorType;


/**
* Function that converts a tensor type to a string.
 */
static const string tensorType2String(const tensorType & type);

/**
* Function that converts a string to a tensor type.
 */
static const tensorType string2TensorType(const string & tensorTypeStr);


private:
		
	/**
	 * The number of columns in the tensor.
	 */
	int NbrOfColumns;
	
	/**
	 * The number of heights of the tensor.
	 */
	int NbrOfHeights;
	
	/**
	 * The number of rows in the tensor.
	 */
	int NbrOfRows;
	
	/**
	 * The number of Widths of the tensor.
	 */
	int NbrOfWidths;
	
	/**
	 * The Type of the tensor.
	 */
	tensorType Type;
	

protected:

	/**
	 * The copy Method.
	 */
	virtual void copy(const tensor & m);

	/**
	 * Method SET for the number of columns.
	 */
	void setNbrOfColumns(const int & nbrOfColumns);
	
	/**
	 * Method SET for the number of heights.
	 */
	void setNbrOfHeights(const int & nbrOfHeights);
	
	/**
	 * Method SET for the number of rows.
	 */
	void setNbrOfRows(const int & nbrOfRows);
	
	/**
	 * Method SET for the number of widths.
	 */
	void setNbrOfWidths(const int & nbrOfWidths);
		
	/**
	 * Method SET for the type of the tensor.
	 */
	void setTensorType(const tensorType & type);
	
public:
	
	/**
	 * The default constructor.
	 */
	tensor(void);
	
	/**
	 * The destructor.
	 */
	virtual ~tensor(void);

	/**
	 * Method Clear : this method makes the current object an empty one.
	 */
	void clear(void);
	
	/**
	 * Method to know if the tensor object is empty.
	 */
	bool empty(void) const;
	
	/**
	 * The GET methd for the coefficients of the tensor.
	 *
	 * This is a pure virtual method.
	 */
	virtual const double & getCoefficient(const int & row , const int & column , const int & width , const int & height)  const=0;
	
	/**
	 * Method GET for the number of columns.
	 */
	const int & getNbrOfColumns(void) const;
	
	/**
	 * Method GET for the number of heights.
	 */
	const int & getNbrOfHeights(void) const;
	
	/**
	 * Method GET for the number of rows.
	 */
	const int & getNbrOfRows(void) const;
	
	/**
	 * Method GET for the number of heights.
	 */
	const int & getNbrOfWidths(void) const;
	
	/**
	 * Method GET for the type of tensor.
	 */
	const tensorType & getTensorType(void) const;
};

/**
* The operator to pretty print the tensor types.
 */
extern ostream & operator<<(ostream & outStream , const tensor::tensorType & type);


#endif
