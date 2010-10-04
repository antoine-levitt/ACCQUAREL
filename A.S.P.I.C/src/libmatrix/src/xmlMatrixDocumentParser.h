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
#ifndef _XML_MATRIX_DOCUMENT_PARSER_
#define _XML_MATRIX_DOCUMENT_PARSER_

#include <xmlDocumentParser.h>
#include "xmlMatrixParser.h"

class xmlMatrixDocumentParser : public xmlDocumentParser , public xmlMatrixParser
{
private:

protected:

	/**
	 * Method GET for the location of the schema.
	 */
	virtual string getSchemaURI(void) const;
	
public:
	
	/**
	 * Default constructor.
	 */
	 xmlMatrixDocumentParser(void);

	/**
	 * Destructor.
	 */
	virtual ~xmlMatrixDocumentParser(void);

	/**
	 * Method to read a matrix symletric in the document.
	 */
	matrixSymetric getMatrixSymetricFull4Document(const string & xmlFileName);

	/**
	 * Method to read a matrix symletric in the document.
	 */
	matrixFull getMatrixFull4Document(const string & xmlFileName);

	/**
	 * The GET Method.
	 */
	//matrix * getMatrix(const string & fileName);
};


#endif
