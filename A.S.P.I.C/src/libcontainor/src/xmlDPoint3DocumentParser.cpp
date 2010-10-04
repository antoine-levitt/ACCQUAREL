/* 
* The containor library of A.S.P.I.C. 
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
#include "xmlDPoint3DocumentParser.h"
#include <aspicConfiguration.h>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlDPoint3DocumentParser::xmlDPoint3DocumentParser(void)
	: xmlDocumentParser()
{
	;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Destructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlDPoint3DocumentParser::~xmlDPoint3DocumentParser(void)
{
	;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the schema URI.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string xmlDPoint3DocumentParser::getSchemaURI(void) const
{
	string schemaURI =  aspicConfiguration::getAspicRoot();
	schemaURI += "/libxmlParser/data/dpoint3.xsd";
	return schemaURI;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET value of the dpoint<3> contained in the document.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dpoint<3> xmlDPoint3DocumentParser::getDPoint34Document(const string & xmlDPoint3FileName)
{
	load(xmlDPoint3FileName);
	return xmlDPoint3Parser::getDPoint3();
}

