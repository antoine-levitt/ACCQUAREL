/* 
* The chemics library of A.S.P.I.C. 
 * Written and directed by Franois Lodier francois.lodier@gmail.com.
 *
 * Copyright (C) 2005  Franois Lodier
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
#include "xmlContractionsParser.h"

//////////////////////////////////////////////////////////////////////
// Chaine de charactres qui contient le nom pour les composantes de
// base des contractions.
///////////////////////////////////////////////////////////////////////
const string xmlContractionsParser::getContractionsListTagName(void)
{
	return "ContractionsList";
}

//////////////////////////////////////////////////////////////////////
// Chaine de charactres qui contient le nom pour les composantes de
// base des contractions.
///////////////////////////////////////////////////////////////////////
const string xmlContractionsParser::getContractionTagName(void)
{
	return "Contraction";
}

////////////////////////////////////////////////////////////////////////
// Chaine de characteres qui contient le nom pour les Exposants de la
// contractions.
/////////////////////////////////////////////////////////////////////////
const string xmlContractionsParser::getExponentTagName(void)
{
	return "Exponent";
}

//////////////////////////////////////////////////////////////////////////
// Chaine de charactres qui contient le nom pour les coefficients de la 
// contractions.
///////////////////////////////////////////////////////////////////////////
const string xmlContractionsParser::getCoefficientTagName(void)
{
	return "Coefficient";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructeur.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlContractionsParser::xmlContractionsParser(void)
: xmlParser()
{
	;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructeur.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlContractionsParser::xmlContractionsParser(DOMElement * xmlContrationElement)
: xmlParser()
{
	setRootElement(xmlContrationElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Destructeur.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlContractionsParser::~xmlContractionsParser(void)
{
	;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for all the contractions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
contractions xmlContractionsParser::getContractions(const DOMElement * rootElement) const
{

	contractions c;
	int contraction , nbrOfContractions;
	DOMElement * contractionsListElement;
	DOMNodeList * contractionNodesList;
	DOMNode * exponentNode;
	DOMNode * coefficientNode;

	/////////////////////////////////////////////////////////////////////////////
	// Rcupration de la racine qui contient les contractions.
	// Lorsque cet lment n'est pas trouv alors on arrete le programme.
	/////////////////////////////////////////////////////////////////////////////
	contractionsListElement = (DOMElement *) getElementByTagName(rootElement  , getContractionsListTagName());

	if(contractionsListElement == NULL) {
		cerr << "Error : in contractions xmlContractionsParser::getContractions(const DOMElement * rootElement) const" << endl;
		cerr << "Error : no contraction list was found in the document." << endl;
		cerr << "Error : aborting." << endl;
		exit(1);	
	}

	//////////////////////////////////////////////////////////////////////////////////
	// Rcupration de la liste des noeuds qui contiennent une contraction.
	//
	// Lorsque cette liste nexiste pas ou est vide on arrete le programme.
	///////////////////////////////////////////////////////////////////////////////////
	contractionNodesList = getElementsByTagName(contractionsListElement,getContractionTagName());

	if(contractionNodesList == NULL || (nbrOfContractions=contractionNodesList->getLength()) == 0) {
		cerr << "Error : in contractions xmlContractionsParser::getContractions(const DOMElement * rootElement) const" << endl;
		cerr << "Error : no contraction  was found in contraction list." << endl;
		cerr << "Error : aborting." << endl;	
		exit(1);
	}

	c.setNbrOfContractions(nbrOfContractions);

	//////////////////////////////////////////////////////////////////////////////////////
	// Ensuite pour chaque contraction on rcupere le coefficient et l'exposant.
	//
	// Evidemment s'il manque l'un des deux on arette le programme.
	//////////////////////////////////////////////////////////////////////////////////////
	for(contraction = 0 ; contraction < nbrOfContractions ; contraction++) {
		
		////////////////////////////////////////////////////////////////////////////////////
		// Rcupration de l'exposant.
		////////////////////////////////////////////////////////////////////////////////////
		exponentNode = getElementByTagName((DOMElement *)contractionNodesList->item(contraction) , getExponentTagName());

		if(exponentNode == NULL) {
			cerr << "Error : in contractions xmlContractionsParser::getContractions(const DOMElement * rootElement) const" << endl;
			cerr << "Error : no exponent  was found for contraction " << contraction << "." << endl;
			cerr << "Error : aborting" << endl;
			exit(1);
		}

		c.setExponent(contraction,getNodeDoubleValue(exponentNode));

		////////////////////////////////////////////////////////////////////////////////////
		// Rcupration du coefficient.
		////////////////////////////////////////////////////////////////////////////////////
		coefficientNode = getElementByTagName((DOMElement *)contractionNodesList->item(contraction) , getCoefficientTagName());

		if(coefficientNode == NULL) {
			cerr << "Error : in contractions xmlContractionsParser::getContractions(const DOMElement * rootElement) const" << endl;
			cerr << "Error : no ceofficient  was found for contraction " << contraction << "." << endl;
			cerr << "Error : aborting" << endl;
			exit(1);
		}

		c.setCoefficient(contraction,getNodeDoubleValue(coefficientNode));
	}

	///////////////////////////////////////////////////////////////////
	// On renvoie les contractions.
	///////////////////////////////////////////////////////////////////
	return c;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for all the contractions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
contractions xmlContractionsParser::getContractions(void) const
{
	return getContractions(getRootElement());
}
