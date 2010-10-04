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
#ifndef _XML_MATRIX_PARSER_
#define _XML_MATRIX_PARSER_

#include "matrix.h"
#include "matrixSymetric.h"
#include "matrixFull.h"

#include <xmlParser.h>

/**
 * This is a very short class to store coeffcient information in
 * an object.
 */
class matrixCoefficientData 
{

public:
	
	/**
	 * the column of the matrix coeffcient.
	 */
	int Column;
	
	/**
	 * The row of the matrix coefficient.
	 */
	int Row;
	
	/**
	 * The value of the matrix coefficient.
	 */
	double Value;
	
};


/**
 * This is the parser of an XML document containing a matrix.
 */
class xmlMatrixParser : virtual public xmlParser
{

private :

protected :

	/**
	 * Method that builds the matrix attempting to 
	 * be a full symetric matrix.
	 */
	//matrix * getMatrixSymetricFull(void) const;

	
public :

	/**
	 * The default contructor.
	 */
	xmlMatrixParser(void);
	
	/**
	 * The default contructor.
	 */
	xmlMatrixParser(DOMElement * rootElement);
	
	/**
	 * the destructor.
	 */
	virtual ~xmlMatrixParser(void);

	/**
	 * Method GET for the tag name containing the coefficient.
	 */
	matrixCoefficientData getCoefficient(int item) const;
	
	/**
	 * Method GET for the tag name containing the coefficient.
	 */
	double getCoefficientsDefaultValue(void) const;
	
	/**
	 * Method GET for the tag name containing the coefficients default value.
	 */
	static const string getCoefficientsDefaultValueTagName(void);
	
	/**
	 * Method GET for the tag name containing the coefficient.
	 */
	static const string getCoefficientTagName(void);
	
	/**
	 * Method GET for the tag name containing the column of a coefficient.
	 */
	static const string getCoefficientColumnTagName(void);	

		/**
		* Method GET for the tag name containing the column of a coefficient.
		 */
	static const string getCoefficientRowTagName(void);	

		/**
		* Method GET for the tag name containing the column of a coefficient.
		 */
	static const string getCoefficientValueTagName(void);	
		
	/**
	 * Method GET for a matrix Symetric with full storage.
	 */
	matrixSymetric getMatrixSymetricFull(void) const;

	/**
	 * Method that builds the matrix attempting to 
	 * be a full symetric matrix.
	 */
	matrixFull getMatrixFull(void) const;

	/**
	 * Method GET for the matrix type.
	 */
	matrix::matrixType getMatrixType(void) const;
	
	/**
		* Method GET for the matrix type.
	 */
	string getMatrixTypeString(void) const;
	
	/**
		* Method to get the tag name for the strucure information.
	 */
	static const string getMatrixTypeAttributeName(void);
	
	/**
	 * Method GET for the numbers of coefficients that are in the file.
	 */
	int getNbrOfCoefficients(void) const;
	
	/**
	 * Method GET for the number of columns.
	 */
	int getNbrOfColumns(void) const;	
	
	/**
	 * Method GET for the tag name containing the number of columns.
	 */
	static const string getNbrOfColumnsTagName(void);	
	
	/**
	 * Method GET for the tag name containing the number of columns.
	 */
	int getNbrOfRows(void) const;	

	/**
	 * Method GET for the tag name containing the number of columns.
	 */
	static const string getNbrOfRowsTagName(void);		
};

#endif
