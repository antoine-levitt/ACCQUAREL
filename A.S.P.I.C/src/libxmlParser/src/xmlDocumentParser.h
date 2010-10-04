/* 
* The xml parser library of A.S.P.I.C. 
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
#ifndef _XML_DOCUMENT_PARSER_
#define _XML_DOCUMENT_PARSER_

#include "xmlParser.h"


/**
 * CLasse générique pour lire les documents XML (DOMDocument).
 *
 * Cette classe permet de lire un fichier XML (méthode load) et de vérifier
 * si cette méthode est compatible avec le schéma qui lui est associé (getSchemaURI).
 */
class xmlDocumentParser : virtual public xmlParser
{

private:
		
	/**
	 * Le parser XERCES qui construit le document DOM associé au document.
	 */
	XercesDOMParser * Parser; 

	/**
	 * Le DOM Document associé au fichier XML.
	 */
	DOMDocument * Document;

protected:
	/**
	 * Method GET for the DOM document.
	 */
	const DOMDocument * getDocument(void) const;

	/**
	 * Method to find the schema associated with the XML
	 * file to parse.
	 */
	virtual string getSchemaURI(void) const =0;
	
	/**
	 * Method SET for the schema.
	 *
	 * The XML Schema Recommendation explicitly 
	 * states that the inclusion of schemaLocation/ 
	 * noNamespaceSchemaLocation attributes in the instance document 
	 * is only a hint; it does not mandate that these attributes must 
	 * be used to locate schemas. This property allows the user to specify 
	 * the no target namespace XML Schema Location externally. If 
	 * specified, the instance document's noNamespaceSchemaLocation 
	 * attribute will be effectively ignored. 
	 *
	 * @param schemaURI The syntax is the same as for the noNamespaceSchemaLocation 
	 * attribute that may occur in an instance document: e.g."file_name.xsd"
		*/
	void setExternalNoNamespaceSchemaLocation(const string & schemaFileName);
		
public :
	
	/**
	 * Constructor.
	 */
	xmlDocumentParser(void);

	/**
	 * Destructor.
	 */
	virtual ~xmlDocumentParser(void);

	/**
	 * Method that parse the document and create the DOM document.
	 */
	void load(const string & xmlFileURI);

	/**
	 * Method that close (destroy) the document.
	 */
	void close(void);
};

#endif


