/* 
* The chemics library of A.S.P.I.C. 
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

#include "xmlBasisElementParser.h"
#include "xmlContractionsParser.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructeur.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlBasisElementParser::xmlBasisElementParser(void)
: xmlParser()
{
	;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructeur avec spécification de la racine.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlBasisElementParser::xmlBasisElementParser(DOMElement * xmlBasisElement)
: xmlParser(xmlBasisElement)
{
	;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Destructeur.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlBasisElementParser::~xmlBasisElementParser(void)
{
	;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is the tag name that encapsulates the basis data for a specific 
// basis element.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlBasisElementParser::getBasisElementTagName(void)
{
	return "BasisElement";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is the name of the Id Attribute to search in the data base with key.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlBasisElementParser::getElementKeyAttributeName(void)
{
	return "ElementKey";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is the name of the Id Attribute to search in the data base with key.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlBasisElementParser::getBasisKeyAttributeName(void)
{
	return "BasisKey";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is the name of the Id Attribute to search in the data base with key.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlBasisElementParser::getShellTypeKeyTagName(void) 
{
	return "ShellTypeKey";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is the name of the Id Attribute to search in the data base with key.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlBasisElementParser::getShellTagName(void) 
{
	return "Shell";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is the name of the Id Attribute to search in the data base with key.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlBasisElementParser::getShellsListTagName(void) 
{
	return "ShellsList";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for a gaussian basis element.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
basisElement xmlBasisElementParser::getBasisElement(const DOMElement * rootElement) const
{
	basisElement element;	
	
	int shell , nbrOfShells;
	DOMElement * basisElementElement;
	DOMElement * shellsListElement;
	DOMNodeList * shellNodesList;
	DOMNode * shellTypeKeyNode;
	xmlContractionsParser contractionsParser;

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// On récupère l'élément de base ...
	//
	// Si on ne trouve pas le tag on arrete le programme.
	// Si on le trouve on attrape la clé de l'élément qui l'identifie. Lorsque l'attribut n'est pas défini
	// cela n'est pas très grave ici, donc on laisse couler.
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	basisElementElement = (DOMElement *) getElementByTagName(rootElement , getBasisElementTagName());

	if(basisElementElement == NULL) {
		cerr << "Error : in basisElement xmlBasisElementParser::getBasisElement(const DOMElement * rootElement) const" << endl;
		cerr << "Error : no basis element was found in the document." << endl;
		cerr << "Error : aborting." << endl;
		assert(0);
		exit(1);
	}

	element.setBasisElementKey(getAttributeValue(basisElementElement , getElementKeyAttributeName()));
	element.setBasisName(getAttributeValue(basisElementElement , getBasisKeyAttributeName()));
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Un élément de base est une liste de couche on cherche donc la liste de couche.
	// 
	// Lorsque cette liste de couche n'est pas trouvée, on arrete le programme.
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	shellsListElement = (DOMElement * ) getElementByTagName(basisElementElement , getShellsListTagName());

	if(shellsListElement == NULL) {
		cerr << "Error : in basisElement xmlBasisElementParser::getBasisElement(const DOMElement * rootElement) const" << endl;
		cerr << "Error : no shells list was found in the basis element." << endl;
		cerr << "Error : aborting." << endl;
		assert(0);
		exit(1);
	}
		
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Ensuite on retrouve la liste des shells.
	//
	// Lorsque cette liste ne compostre aucun élément on arrete le programme.
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	shellNodesList = getElementsByTagName(shellsListElement , getShellTagName());

	if(shellNodesList == NULL || (nbrOfShells = shellNodesList->getLength()) == 0) {
		cerr << "Error : in basisElement xmlBasisElementParser::getBasisElement(const DOMElement * rootElement) const" << endl;
		cerr << "Error : no shells was found in the shells list." << endl;
		cerr << "Error : aborting." << endl;
		assert(0);
		exit(1);
	}
	
	element.setNbrOfShells(nbrOfShells);

	///////////////////////////////////////////////////////////////////////////////////////////
	// Pour chaque couche on extrait ensuite les informations nécessaires 
	// à sa construction.
	////////////////////////////////////////////////////////////////////////////////////////////
	for(shell = 0; shell < nbrOfShells ; shell++) {
		
		//////////////////////////////////////////////////////////////////////////////
		// On récupère le noued qui contient le type de couche.
		//
		// -On arrete le programme lorsque l'on  ne trouve pas le type de
		// couche mais aussi lorsque l'ontrouve un type qui n'est pas valide.
		//////////////////////////////////////////////////////////////////////////////
		shellTypeKeyNode = getElementByTagName((DOMElement *) shellNodesList->item(shell),getShellTypeKeyTagName());
		
		if(shellTypeKeyNode == NULL) {
			cerr << "Error : in basisElement xmlBasisElementParser::getBasisElement(const DOMElement * rootElement) const" << endl;
			cerr << "Error : no shell type  was found in the shell " << shell << " in the list." << endl;
			cerr << "Error : aborting." << endl;
			assert(0);
			exit(1);
		}

		element.setShellType(shell , getNodeStringValue(shellTypeKeyNode,Remove_White_Space));
		
		if(element.getShellType(shell).empty()) {
			cerr << "Error : in basisElement xmlBasisElementParser::getBasisElement(const DOMElement * rootElement) const" << endl;
			cerr << "Error : no shell type  was defined for key \"" << getNodeStringValue(shellTypeKeyNode,Remove_White_Space) << "\"." << endl;
			cerr << "Error : aborting." << endl;
			assert(0);
			exit(1);		
		}

		//////////////////////////////////////////////////////////////////////////////
		// Dans le cas des contractions on refourgue le bébé au parser 
		// approprié.
		///////////////////////////////////////////////////////////////////////////////
		element.setContractions4Shell(shell , contractionsParser.getContractions((DOMElement *) shellNodesList->item(shell)));
	
		if(element.getContractions4Shell(shell).empty()) {
			cerr << "Error : in basisElement xmlBasisElementParser::getBasisElement(const DOMElement * rootElement) const" << endl;
			cerr << "Error : no contractions were defined for key \"" << getNodeStringValue(shellTypeKeyNode,Remove_White_Space) << "\"." << endl;
			cerr << "Error : aborting." << endl;
			assert(0);
			exit(1);		
		}

	}


	/////////////////////////////////////////////////////////
	// On renvoie l'élément de base ainsi construit.
	/////////////////////////////////////////////////////////
	return element;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for a gaussian basis element.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
basisElement xmlBasisElementParser::getBasisElement(void) const
{
	return getBasisElement(getRootElement());
}
