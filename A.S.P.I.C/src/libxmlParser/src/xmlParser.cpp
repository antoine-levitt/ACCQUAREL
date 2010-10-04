/* 
* The xml parser library of A.S.P.I.C. 
 * Written and directed by Fran�ois Lodier francois.lodier@gmail.com.
 *
 * Copyright (C) 2005  Fran�ois Lodier
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
#include "xmlParser.h"


////////////////////////////////////////////////////////////////////////////////////////////
// Construcuteur.
////////////////////////////////////////////////////////////////////////////////////////////
xmlParser::xmlParser(void)
: RootElement(NULL)
{
	;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Construcuteur.
////////////////////////////////////////////////////////////////////////////////////////////
xmlParser::xmlParser(DOMElement * rootElement)
: RootElement(NULL)
{
	setRootElement(rootElement);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Destructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlParser::~xmlParser(void)
{
	RootElement = NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for an attribute value.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string xmlParser::getAttributeValue(DOMElement * element , const string & attributeName) const
{

	XMLCh * xmlAttributeName = XMLString::transcode(attributeName.c_str());
	const XMLCh * xmlAttributeValue = element->getAttribute(xmlAttributeName);
	XMLString::release(&xmlAttributeName);

	if(xmlAttributeValue == NULL) {
		return "";
	}
	
	char * charAttributeValue = XMLString::transcode(xmlAttributeValue);
	string strAttributeValue = charAttributeValue;
	
	XMLString::release(&charAttributeValue);
	
	return strAttributeValue;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for an attribute value.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string xmlParser::getAttributeValue(const string & attributeName) const
{
	
	if(getRootElement() == NULL) {
			//cerr << "Warning : in DOMNode * xmlParser::getElementByTagName(const string & elementName) const" << endl;
			//cerr << "Warning : searching for tags in an empty element" << endl;
			//cerr << "Warning : a Null pointer will be return" << endl;
			return "";
	}
	
	return getAttributeValue(getRootElement() , attributeName);
}

////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour accéder à un élément unique dont on connait le nom 
// dans le fichier XML.
//
// Description : Cette méthode permet de retrouver la partie <elementName> ... </elementName>
// du duocment XML.
//
// Note : Si l'élément elementName n'est pas trouvé, on renvoie le pointeur nul.
//
// Note : Si plusieurs éléments sont trouvés, on renvoie le premier.
/////////////////////////////////////////////////////////////////////////////////////////////
DOMNode * xmlParser::getElementByTagName(const DOMElement * xmlElement , const string & elementName) const
{
	DOMNodeList * elementList = getElementsByTagName(xmlElement , elementName);

	if(elementList == NULL || elementList->getLength() == 0) {
		//cerr << "Error : in DOMNode * gaussianBasisFunctionParser::getElementByTagName(const string & elementName) const" << endl;
		//cerr << "Error : no element with tag name " << elementName << " where found." << endl;
		//cerr << "Error : returning null" << endl;
		
		char * chLocalName = XMLString::transcode(xmlElement->getLocalName());
		string localName;
		if(chLocalName == NULL) {
			localName = "";
		} else {
			localName = chLocalName;
		}
		XMLString::release(&chLocalName);
	
		if(localName == elementName) 
			return (DOMNode *)xmlElement;
		else 
			return NULL;
	}
	
	if(elementList->getLength() > 1) {
		cerr << "Warning : in DOMNode * gaussianBasisFunctionParser::getElementByTagName(const string & elementName) const" << endl;
		cerr << "Warning : element with tag name " << elementName << " is not unique" << endl;
		cerr << "Warning : returning first element of the list." << endl;
	}

	return elementList->item(0);	
}


//////////////////////////////////////////////////////////////////////////
// Méthode pour accéder au noeud elementName dans cette partie de 
// l'arboréscence.
///////////////////////////////////////////////////////////////////////////
DOMNode * xmlParser::getElementByTagName(const string & elementName) const
{
	if(getRootElement() == NULL) {
		//cerr << "Warning : in DOMNode * xmlParser::getElementByTagName(const string & elementName) const" << endl;
		//cerr << "Warning : searching for tags in an empty element" << endl;
		//cerr << "Warning : a Null pointer will be return" << endl;
		return NULL;
	}

	return getElementByTagName(getRootElement(),elementName);
}

////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour accéder à une liste d'éléments dont on connait le nom 
// dans le fichier XML.
//
// Description : Cette méthode permet de retrouver toutes les parties <elementName> ... </elementName>
// du duocment XML.
//
// Note : Si l'élément elementName n'est pas trouvé, on renvoie le pointeur nul.
//
// Note : Si plusieurs éléments sont trouvés, on renvoie le premier.
/////////////////////////////////////////////////////////////////////////////////////////////
DOMNodeList * xmlParser::getElementsByTagName(const DOMElement * xmlElement , const string & elementName) const 
{
	if(xmlElement == NULL) {
		//cerr << "Warning : in DOMNodeList * xmlParser::getElementsByTagName(const DOMElement * xmlElement , const string & elementName) const" << endl;
		//cerr << "Warning : searching for tags in an empty element" << endl;
		//cerr << "Warning : a Null pointer will be return" << endl;
		return NULL;
	}

	XMLCh * xmlElementName = XMLString::transcode(elementName.c_str());
	DOMNodeList * elementList = xmlElement->getElementsByTagNameNS(xmlElement->getNamespaceURI(),xmlElementName);
	XMLString::release(&xmlElementName);
	return elementList;	
}

///////////////////////////////////////////////////////////////////////////////
// Méthode pour chercher un élément dans la partie de l'arboréscence
// que l'on manipule avec le xmlParser.
////////////////////////////////////////////////////////////////////////////////
DOMNodeList * xmlParser::getElementsByTagName(const string & elementName) const 
{
	if(RootElement == NULL) {
		//cerr << "Warning : in DOMNodeList * xmlParser::getElementsByTagName(const string & elementName) const" << endl;
		//cerr << "Warning : searching for tags in an empty element" << endl;
		//cerr << "Warning : a Null pointer will be return" << endl;
		return NULL;
	}

	return getElementsByTagName(RootElement , elementName);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// M�thode pour connaitre le nombre de fils d'un �lement xmlElement dont le tag est elementName.
////////////////////////////////////////////////////////////////////////////////////////////////////////////
int xmlParser::getNbrOfElementsWithTagName(const DOMElement * xmlElement , const string & elementName) const
{
	DOMNodeList * elementList = getElementsByTagName(xmlElement , elementName);
	if(elementList == NULL) {
		return 0;
	}	
	return elementList->getLength();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// M�thode pour connaitre le nombre de fils d'un �lement xmlElement dont le tag est elementName.
////////////////////////////////////////////////////////////////////////////////////////////////////////////
int xmlParser::getNbrOfElementsWithTagName(const string & elementName) const
{
	
	if(RootElement == NULL) {
		return 0;
	}
	
	return getNbrOfElementsWithTagName(RootElement , elementName);
}

////////////////////////////////////////////////////////////////
// Méthode qui extrait la valeur double d'un neud.
////////////////////////////////////////////////////////////////
double xmlParser::getNodeDoubleValue(const DOMNode * node) const
{
	if(node == NULL) {
		//cerr << "Warning : in double xmlParser::getNodeDoubleValue(const DOMNode * node) const" << endl;
		//cerr << "Warning : cannot acces the node value of an empty node." << endl;
		//cerr << "Warning : returning NaN." << endl;
		return 0; 
	}

	XMLDouble xmlDoubleValue(getNodeTextValue(node));
	return xmlDoubleValue.getValue();
}

//////////////////////////////////////////////////////////////
// M�thode pour connaitre la valeur d'un noued sous forme
// d'entier.
//
// @warning lorsque le pointeur pass� en argument est nul
// alors la valeur retourn�e est nan.
//////////////////////////////////////////////////////////////
int xmlParser::getNodeIntegerValue(const DOMNode * node) const
{
	if(node == NULL) {
		//cerr << "Warning : in double xmlParser::getNodeIntegerValue(const DOMNode * node) const" << endl;
		//cerr << "Warning : cannot acces the node value of an empty node." << endl;
		//cerr << "Warning : returning NaN." << endl;
		return 0; 
	}

	XMLBigInteger xmlInteger(getNodeTextValue(node));
	return xmlInteger.intValue();
}

/////////////////////////////////////////////////////////////////////
// M�thode pour connaitre la valeur d'un �lement sous forme d'une
// chaine de charct�re.
//
// @warning Lorsque la valeur pass�e en argument est NULL alors 
// on affiche un message d'erreur et on renvoie la chaine vide.
/////////////////////////////////////////////////////////////////////
string xmlParser::getNodeStringValue(const DOMNode * node , stringFormatOption opt) const
{
	string value;
	
	if(node == NULL) {
		//cerr << "Error : in double xmlParser::getNodeStringValue(const DOMNode * node) const" << endl;
		//cerr << "Error : the node argument is NULL." << endl;
		//cerr << "Error : no value ca be found. Return 0." << endl;
		value = "";
		return value;
	}
	
	
	char * charValue;
	XMLCh * xmlValue = XMLString::replicate( getNodeTextValue(node) );
	
	// Suivant l'option de formattage on va changer
	// ce qu'il faut dans la chaine de charact�re.
	switch(opt) {
		
		case No_Format :
			charValue = XMLString::transcode(xmlValue);
			break;
			
		case Remove_White_Space:
			XMLString::removeWS(xmlValue);
			charValue = XMLString::transcode(xmlValue);
			break;
			
		case Collapse_White_Space : 
			XMLString::collapseWS(xmlValue);
			charValue = XMLString::transcode(xmlValue);
			break;
			
		default:
			cerr << "Error : in " << endl;
			cerr << "Error : the option to format the string was not recognise." << endl;
			cerr << "Error : nothing will be done" << endl;
			return "";
	}
	
	// On trasforme ce qu'il faut pour que cela fonctionne.
	value = charValue;
	
	// On nettoye la m�moire allou�e.
	XMLString::release(&charValue);
	XMLString::release(&xmlValue);
	// Et on renvoie la valeur.
	return value;
}

/////////////////////////////////////////////////////////////////////
// Méthode pour connaitre la vaeur d'un element sous forme de
// chaine de charactère.
/////////////////////////////////////////////////////////////////////
const XMLCh * xmlParser::getNodeTextValue(const DOMNode * node) const
{
	if(node == NULL) {
		//cerr << "Error : in double xmlParser::getNodeTextValue(const DOMNode * node) const" << endl;
		//cerr << "Error : the node argument is NULL." << endl;
		//cerr << "Error : no value ca be found. Return 0." << endl;
		return 0;
	}

	DOMNode * nodeValue = node->getFirstChild();
	
	if(nodeValue == NULL || nodeValue->getNodeType() != DOMNode::TEXT_NODE) { 
		//cerr << "Error : in double xmlParser::getNodeDoubleValue(const DOMNode * node) const" << endl;
		//cerr << "Error : the node that you pass as an argument does not seem to contain a double" << endl;
		//cerr << "Error : returning 0." << endl;
		return 0;
	}

	return nodeValue->getNodeValue();
}

/////////////////////////////////////////////////////////////////////////
// Method GET for the root element.
/////////////////////////////////////////////////////////////////////////
DOMElement * xmlParser::getRootElement(void) const
{
	return RootElement;
}
/////////////////////////////////////////////////////////////////////////
// Method SET for the root element.
/////////////////////////////////////////////////////////////////////////
void xmlParser::setRootElement(DOMElement * rootElement)
{
	if(rootElement == NULL) {
		//cerr << "Warning : in void xmlParser::setRootElement(const Element * rootElement)" << endl;
		//cerr << "Warning : setting empty root element" << endl;
		//cerr << "Warning : no data will be accessible." << endl;
	}

	RootElement = rootElement;
}

