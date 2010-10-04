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
#ifndef _TENSOR_SYMETRIC_
#define _TENSOR_SYMETRIC_


#include <fstream>
#include <iomanip>
#include <iostream>
#include "tensorNonSparse.h"

using namespace std;

class tensorSymetric : public tensorNonSparse
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
	virtual int hash( const int & row , const int & column , const int & width , const int & height) const;
	
public:
		
	/**
		* The default constructor.
		*/
		tensorSymetric(void);
	
	/**
		* The constructor with size specification.
	 *
	 * This method constructs a symetric tensor that will be aible to
	 * store nbrOfRows x nbrOfRows elements.
	 */
	tensorSymetric(const int & nbrOfRows);
	
	/**
		* The destructor.
	 */
	~tensorSymetric(void);	
	
	/**
		* Method read : This method reads a tensor in
	 * an unformated file.
	 */
	void readText(const string & fileName);
	
	/**
		* Method read : This method reads a tensor in
	 * an unformated file.
	 */
	void readText(istream & inStream);
	
	/**
		* Method SET for the tensor Size.
	 */
	void setTensorSize(int nbrOfRows);
	
	/**
		* Method SET for the tensor Size.
	 *
	 * @warning this will fail if the numbers of rows
	 * is different from the number of columns.
	 */
	void setTensorSize(const int & nbrOfRows , const int & nbrOfColumns , const int & nbrOfWidth , const int & nbrOfHeights);
	
	/**
		* Method that writes a tensor symetric in a stream.
	 */
	void writePretty(ostream & outStream , int width = 10 , int precision = 6) const;
	
	/**
		* Method that writes the symetric Tensor in a file for that can be 
	 * read with the method readText.
	 */
	void writeText(ostream & outStream , int precision = 6) const;

};

#endif
