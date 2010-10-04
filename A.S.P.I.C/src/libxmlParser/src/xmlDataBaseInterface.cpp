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
#include <assert.h>
#include "xmlDataBaseInterface.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Le constructeur
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlDataBaseInterface::xmlDataBaseInterface(const string & entryTagName , const string & idAttributeName)
:	EntryTagName(entryTagName) , IdAttributeName(idAttributeName)
{
	;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Le destructeur.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlDataBaseInterface::~xmlDataBaseInterface(void)
{
	;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode GET pour une entée dans la base de donnée
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
DOMNode * xmlDataBaseInterface::getDataBaseEntry(int item) const
{
	assert(item >=0);
	assert(item < getNbrOfEntries());
	
	return getElementsByTagName(EntryTagName)->item(item);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the ID of an entry of the data base.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
string xmlDataBaseInterface::getID4Entry(int item) const
{
	assert(item >= 0);
	assert(item < getNbrOfEntries());
	
	DOMElement * entryElement =(DOMElement *) getDataBaseEntry(item);
	XMLCh * xmlIdAttributeName;
	char * charId;
	string strId;

	xmlIdAttributeName = XMLString::transcode(IdAttributeName.c_str());
	charId = XMLString::transcode(entryElement->getAttribute(xmlIdAttributeName));
	XMLString::release(&xmlIdAttributeName);
	strId = charId;
	XMLString::release(&charId);
	return strId;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour retrouver les information d'un élément à partir  de sa clé dans la base de données.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
DOMElement * xmlDataBaseInterface::getElementById(const string & key) const
{
	XMLCh * xmlKey = XMLString::transcode(key.c_str());
	XMLString::removeWS(xmlKey);

	if(XMLString::stringLen(xmlKey) == 0) {
		cerr << "Error : in DOMElement * xmlDataBaseInterface::getElementById(const string & key) const" << endl;
		cerr << "Error : the key argument is empty." << endl;
		cerr << "Error : return NULL" << endl;
		return NULL;
	}
	
	DOMElement * element = getDocument()->getElementById(xmlKey);
	XMLString::release(&xmlKey);
	
	if(element == NULL) {
		cerr << "Error : in DOMElement * xmlDataBaseInterface::getElementById(const string & key) const" << endl;
		cerr << "Error : no element with key : \"" << key << "\" was found." << endl;
		cerr << "Error : return NULL" << endl;
	}
	
	return element;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode GET pour le nombre d'éléments que contient la base de données.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int xmlDataBaseInterface::getNbrOfEntries(void) const
{
	return getNbrOfElementsWithTagName(EntryTagName);
}
