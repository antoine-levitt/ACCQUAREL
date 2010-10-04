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

#include "tensor.h"
#include "tensorFull.h"
#include "tensorSymetric.h"
#include <xmlParser.h>

/*
 * Class that stores the only information about a coefficient
 * of the tensor.
 */
class tensorCoefficientData
{
private:

protected:

public :
	
	/**
	 * The column of the tensor element.
	 */
	int Column;

	/**
	 * The height of the tensor element.
	 */
	int Height;
	
	/**
	 * The row of the tensor element.
	 */
	int Row;
	
	/**
	 * The width of the tensor element.
	 */
	int Width;
	
	/**
	 * The value of the tensor element.
	 */
	double Value;
};


/*
 * Class that parse a document to read a tensor.
 */
class xmlTensorParser : public virtual xmlParser
{
	
private:

protected:
	/**
	 * Method that builds the tensor attempting to 
	 * be a full symetric tensor.
	 */
	tensor * getTensorSymetricFull(void) const;
	
	/**
	 * Method that builds the tensor attempting to 
	 * be a full symetric tensor.
	 */
	tensor * getTensorFull(void) const;
	
	
public:
	
	/*
	 * The default constructor.
	 */
	xmlTensorParser(void);
	
	/*
	 * The constructor with the root element.
	 */
	xmlTensorParser(DOMElement * rootElement);
	
	/*
	 * the destructor.
	 */
	virtual ~xmlTensorParser(void);
	
	/*
	 * Method GET for the tag name containing the coefficient.
	 */
	tensorCoefficientData getCoefficient(int item) const;
	
	/*
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
	static const string getCoefficientHeightTagName(void);	
	
	/**
	 * Method GET for the tag name containing the column of a coefficient.
	 */
	static const string getCoefficientRowTagName(void);	
	
	/**
	 * Method GET for the tag name containing the column of a coefficient.
	 */
	static const string getCoefficientValueTagName(void);	
	
	/**
	 * Method GET for the tag name containing the column of a coefficient.
	 */
	static const string getCoefficientWidthTagName(void);	
	
	/*
	 * Method GET for the number of coefficients.
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
	int getNbrOfHeights(void) const;	
	
	/**
	 * Method GET for the tag name containing the number of columns.
	 */
	static const string getNbrOfHeightsTagName(void);		
	
	/**
		* Method GET for the tag name containing the number of columns.
	 */
	int getNbrOfRows(void) const;	
	
	/**
		* Method GET for the tag name containing the number of columns.
	 */
	static const string getNbrOfRowsTagName(void);			

	/**
		* Method GET for the number of columns.
	 */
	int getNbrOfWidths(void) const;	
	
	/**
	 * Method GET for the tag name containing the number of columns.
	 */
	static const string getNbrOfWidthsTagName(void);	
	
	/**
	 * Method GET for the tensor.
	 */
	tensor * getTensor(void) const;
	
	/**
	 * Method GET for the tensor type.
	 */
	tensor::tensorType getTensorType(void) const;
	
	/**
	 * Method GET for the tensor type.
	 */
	string getTensorTypeString(void) const;
	
	/**
	 * Method to get the tag name for the strucure information.
	 */
	static const string getTensorTypeAttributeName(void);
	
};


