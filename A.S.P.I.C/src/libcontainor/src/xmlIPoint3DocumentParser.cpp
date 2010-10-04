/* 
* The containor library of A.S.P.I.C. 
 * Written and directed by François Lodier support.aspic@gmail.com.
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
#include "xmlIPoint3DocumentParser.h"
#include <aspicConfiguration.h>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlIPoint3DocumentParser::xmlIPoint3DocumentParser(void)
	: xmlDocumentParser()
{
	;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Destructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlIPoint3DocumentParser::~xmlIPoint3DocumentParser(void)
{
	;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the schema URI.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string xmlIPoint3DocumentParser::getSchemaURI(void) const
{
	string schemaURI =  aspicConfiguration::getAspicRoot();
	schemaURI += "/data/containor/ipoint3.xsd";
	return schemaURI;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET value of the dpoint<3> contained in the document.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ipoint<3> xmlIPoint3DocumentParser::getIPoint34Document(const string & xmlIPointFileName)
{
	load(xmlIPointFileName);
	return xmlIPoint3Parser::getIPoint3();
}

