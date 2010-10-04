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
#include "xmlAtomParser.h"
#include <xmlDPoint3Parser.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour trouver le nom du tag qui va contenir toutes les informations
// sur l'élement chimique.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlAtomParser::getAtomTagName(void)
{
	return "Atom";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour trouver le nom du tag qui va contenir toutes les informations
// sur l'élement chimique.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlAtomParser::getChemicalElementTagName(void)
{
	return "ChemicalElement";
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour trouver le nom du tag qui contient la clé de l'élement chimique.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlAtomParser::getChemicalElementKeyTagName(void)
{
	return "ChemicalElementKey";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour retrouver le nom du tag qui contient le nombre de neutrons 
// de l'élément chimique.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlAtomParser::getNbrOfNeutronsTagName(void)
{
	return "NbrOfNeutrons";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Le nom de la balise qui contient les informations sur le positionnement de l'atome.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlAtomParser::getAtomCenterTagName(void)
{
	return "AtomCenter";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Le nom de la balise qui contient les information à propos des unitées des atomes.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlAtomParser::getAtomCenterUnitTagName(void)
{
	return "DistanceUnit";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor of the class.
//
// @param rootAtomElement the root node that contains all informations about the atom that is gone be read.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlAtomParser::xmlAtomParser(void)
{
	;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor of the class.
//
// @param rootAtomElement the root node that contains all informations about the atom that is gone be read.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlAtomParser::xmlAtomParser(DOMElement * rootAtomElement)
{
	setRootElement(rootAtomElement);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Destructor of the class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlAtomParser::~xmlAtomParser(void)
{
	;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour construire totalement un atome.
//
// @return l'atome initialisé à partir de données du document XML.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
atom xmlAtomParser::getAtom(void) const
{
	atom a;	
	pair<dpoint<3>,atom::distanceUnit> atomCenter;
	pair<chemicalElement,int> atomChemicalElement;
	
	////////////////////////////////////////////////////////////////////////
	// On commence par l'élément chimique sur lequel
	// s'appuie l'atome.
	// 
	// - Récupération des informations concernant l'élément chimique.
	// - On recopie l'élément chipmique trouvé.
	// - On recopie le nombre de neutrons trouvé.
	////////////////////////////////////////////////////////////////////////
	atomChemicalElement = getAtomChemicalElement();
	a.setChemicalElement(atomChemicalElement.first);
	a.setNbrOfNeutrons(atomChemicalElement.second);

	////////////////////////////////////////////////////////////////////////
	// On Gère ensuite ce qui a trait à la positon de l'atome.
	//
	// - Récupération des information sous forme simple.
	// - On recopie la position.
	// - On recopie l'unité.
	/////////////////////////////////////////////////////////////////////////
	atomCenter = getAtomCenter();
	a.setPosition(atomCenter.first);
	a.setDistanceUnit(atomCenter.second);
	return a;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour récupérer les informations qui concerne l'élément chimique constituant 
// l'atome.
//
// @return the key for the chemical element. 
//
// @warning if the chemical element referecence is not found then an empty string is returned.	 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
pair<chemicalElement,int> xmlAtomParser::getAtomChemicalElement(void) const
{
	DOMElement * chemicalElementElement;
	DOMNode * chemicalElementKeyNode;
	DOMNode * chemicalElementNbrOfNeutronsNode;
	pair<chemicalElement,int> atomChemicalElement;
	string chemicalElementKey;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// On commence par récupérer la partie du document qui correspond à la description 
	// de l'élément chimique. En fait cette partie est comprise entre les balises <ChemicalElement>
	// et </ChemicalElement>.
	// 
	// - Lorsque l'on ne trouve pas cela on va interrompre le programme car nous avons quelque
	// chose qui ne va pas.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	chemicalElementElement = (DOMElement *)getElementByTagName(getChemicalElementTagName());

	if(chemicalElementElement == NULL) {
		cerr << "Error : in pair<chemicalElement,int> xmlAtomParser::getAtomChemicalElement(void) const" << endl;
		cerr << "Error : no chemical element was defined." << endl; 
		cerr << "Error : atoms need a chemical element description to be well defined. Aborting." << endl;
		exit(1);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Ensuite on récupère les données pour l'élément chimique à proprement parler.
	// 
	// - On récupère la clé de l'élément chimique.
	// - A l'aide de la base de données des éléments chimiques on complète les informations
	// qui le concerne.
	// 
	// - Lorsque la cké de l'élément chimique n'est pas définie dans le document XML on 
	// arrète le programme.
	// - Lorsque l'élément chimique n'a pas été trouvé dans la base de données on arrete le 
	// programme.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	chemicalElementKeyNode = getElementByTagName(chemicalElementElement , getChemicalElementKeyTagName());

	if(chemicalElementKeyNode == NULL) {
		cerr << "Error : in pair<chemicalElement,int> xmlAtomParser::getAtomChemicalElement(void) const" << endl;
		cerr << "Error : no chemical element key was defined." << endl; 
		cerr << "Error : atoms need a chemical element key to be well defined. Aborting." << endl;	
		exit(1);
	}

	chemicalElementKey = getNodeStringValue(chemicalElementKeyNode , xmlParser::Remove_White_Space);
	atomChemicalElement.first = xmlChemicalElementsDataBaseInterface::getChemicalElement(chemicalElementKey);

	if(atomChemicalElement.first.empty()) {
		cerr << "Error : in pair<chemicalElement,int> xmlAtomParser::getAtomChemicalElement(void) const" << endl;
		cerr << "Error : no chemical element with key \"" << chemicalElementKey << "\" was found in the data base" << endl; 
		cerr << "Error : please check file:///$ASPICROOT/documentation/chemicalElements.html for more informations." << endl;	
		exit(1);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Enfin on récupère le nombre de neutrons associé à cet élément chimique.
	// 
	// - Lorsque le nombre de neutrons n'est pas défini dans le document XML on renvoie l'isotope 
	// le plus probable.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	chemicalElementNbrOfNeutronsNode = getElementByTagName(chemicalElementElement , getNbrOfNeutronsTagName());

	if(chemicalElementNbrOfNeutronsNode == NULL) {
		atomChemicalElement.second = atomChemicalElement.first.getNbrOfNeutrons4Isotope(0);
	} else {
		atomChemicalElement.second = getNodeIntegerValue(chemicalElementNbrOfNeutronsNode);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Pour finir on renvoie les informations récupérées dans le fichier.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return atomChemicalElement;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour récupérer les infirmations sur le centre de l'atome.
// <AtomCenter>
//	
// <Position>
//	<x> ... </x>
//	<y> ... </y>
//	<z> ... </z>
// </Position>
//
// <DistanceUnit> ... <DistanceUnit>
// </AtomCenter> 
//
// @return les information nécessaire au centrage de l'atome, càd une paire composée des coordonnées
// du centre de l'atome et de l'unité dans laquelle celles ci sont exprimées.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
pair<dpoint<3>  , atom::distanceUnit> xmlAtomParser::getAtomCenter(void) const
{
	DOMElement * atomCenterElement;
	DOMNode * unitNode;
	pair<dpoint<3> ,  atom::distanceUnit> atomCenter;
	xmlDPoint3Parser positionParser;
	

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 1 - On récupère la partie du document qui correspond
	// au centre de l'atome.
	//
	// Lorsque cette partie n'est pas trouvée alors on arrete l'execution 
	// du programme car on défini un atome et pas un élément chimique.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	atomCenterElement =  (DOMElement *)getElementByTagName(getAtomCenterTagName());

	if(atomCenterElement == NULL) {
		cerr << "Error : in dpoint<3> xmlAtomParser::getAtomCenterPosition(void) const" << endl;
		cerr << "Error : no atom center was defined." << endl;
		cerr << "Error : atoms need a center description to be well defined. Aborting." << endl;
		exit(1);
	}
		
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 2 - Ici on récupère la partie coordonnées pour le centre de l'atome.
	//
	// - Comme c'est un dpoint<3> que cela doit contenir on cherche le fils avec la racine d'un 
	// xmlDPoint3Parser et on créé le parser à partir de l'élément retrouvé.
	// !! Si on ne trouve pas se fils cela provoque l'arret du programme.
	//
	// - Ensuite on laisse le parser faire son travail.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	atomCenter.first = positionParser.getDPoint3(atomCenterElement);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 3 - Enfin on récupère l'unité dans laquelle sont exprimées les coordonnées de la position.
	//
	//  - On cherche l'élément qui correspond à l'unité. Comme ceci est optionnel lorsque l'on ne le 
	// trouve pas on renvoie l'unité par défaut qui est l'unité atomique.
	// 
	// - Lorsqu'on le trouve on converti la valeur de son contenu en atom::distanceUnit.
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	unitNode = getElementByTagName(atomCenterElement , getAtomCenterUnitTagName());

	if(unitNode == NULL) {
		atomCenter.second = atom::ATOMIC_UNIT;
	} else {
		atomCenter.second = atom::string2DistanceUnit(getNodeStringValue(unitNode , xmlParser::Remove_White_Space));		
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 4 - Il ne reste plus qu'à renvoyer ce que nous avons trouvé.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return atomCenter;
}
