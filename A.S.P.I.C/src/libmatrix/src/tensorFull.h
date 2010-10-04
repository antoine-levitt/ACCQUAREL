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
#ifndef _TENSOR_FULL_
#define _TENSOR_FULL_

#include <iomanip>
#include <iostream>
#include "tensorNonSparse.h"
using namespace std;


class tensorFull : public tensorNonSparse
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
	virtual int hash(const int & row , const int & column , const int & width , const int & height) const;

public:
	
	/**
	 * The constructor.
	 */
	tensorFull(void);

	/**
	 * The constructor.
	 */
	tensorFull( const int & nbrOfRows , const int & nbrOfColumns , const int & nbrOfWidths , const int & nbrOfHeights);
	
	/**
	 * The destructor.
	 */
	~tensorFull(void);	
	
	/**
	 * The method to write tensores in the prettiest way.
	 */
	void writePretty(ostream & outStream , int width = 10 , int precision = 6) const;
};

#endif