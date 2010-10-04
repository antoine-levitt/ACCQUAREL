/* 
 * The chemics library of A.S.P.I.C. 
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
#ifndef _ATOM_PARSER_
#define _ATOM_PARSER_

#include "atom.h"
#include "chemicalElementsDataBaseInterface.h"
#include <dpoint.h>
#include <iostream>
#include <xmlParser.h>
using namespace std;

/**
 * Class Pour créer dans un programme C++ un objet de type
 * atome qui contient les informations contenus dans un fichier
 * XML.
 *
 * Voici un exemple de ce que sait lire la classe xmlAtomParser :
 * 
 *	<Atom>
 *		
 *		<ChemicalElement>
 *			<ChemicalElementKey> H </chemicalElementKey>
 *			<NbrOfNeutrons> 1 </NbrOfNeutrons>
 *   </ChemicalElement>
 * 
 *		<AtomCenter>
 *			<Position> 
 *				<x> 0 </x> 
 *				<y> 0 </y> 
 *				<z> 0 </z> 
 *			</Position>
 *
 *			<DistanceUnit> AU </DistanceUnit>
 *		</AtomCenter>
 *
 *	</Atom>
 *
 * A partir d'un tel fichier cette classe construit un (objet) atome qui contient les informations :
 *
 * - L'élément chimique possède un clé H:
 *		- son symbol est H.
 *	  - son nom est l'Hydrogène.
 *	  - Il contient 1 proton.
 *		- Il contient 1 neutron.
 *  
 * - La position de cet atome est (0,0,0) et elle est exprimée en unités atomiques.
 */
class xmlAtomParser : public xmlParser
{
private:
			
protected:
	
	/**
	 * Méthode pour retrouver les informations qui concernent le 
	 * centre de l'atome.
	 *
	 * @return le couple formée par la position de l'atome et de l'unité dans laquelle
	 * la postion est exprimée.
	 */
	pair<dpoint<3> , atom::distanceUnit> getAtomCenter(void) const;

	/**
	 * Méthode pour retrouver les information relative à l'élément chimique
	 * que l'on manipule.
	 *
	 * En fait cette méthode construit non seulement l'élément chimique sur
	 * lequel s'appuie l'atome mais le spécifie aussi car, pour un atome, le nombre
	 * de neutrons doit etre connu. C'est pourquoi elle renvoie une paire et pas simplement
	 * un élément chimique.
	 *
	 * @return la paire formé par l'élément chimique sur lequel s'appuie l'atome et le nombre
	 * de neutrons de cet élément.
	 */
	pair<chemicalElement,int> getAtomChemicalElement(void) const;
	
public:
	/**
	 * Méthode qui contient le nom du tag racine qui contient un atome.
	 */
	static const string getAtomTagName(void);

	/**
	 * Méthode pour trouver le nom du tag qui va contenir toutes les informations
	 * sur l'élement chimique.
	 */
	static const string getChemicalElementTagName(void);

	/**
	 * Méthode pour trouver le nom du tag qui contient la clé de l'élement
	 * chimique.
	 */
	static const string getChemicalElementKeyTagName();
	
	/**
	 *	Méthode pour retrouver le nom du tag qui contient le nombre de neutrons 
	 * de l'élément chimique.
	 */
	static const string getNbrOfNeutronsTagName();
	
	/**
	 * The name of the tag that contains the position of the atom.
	 */
	static const string getAtomCenterTagName();	

	/**
	 * The name of the tag that contains the position of the atom.
	 */
	static const string getAtomCenterUnitTagName();	

	/**
	 * Constructeur de la classe.
	 */
	xmlAtomParser(void);
	
	/**
	 * Constructeur de la classe.
	 *
	 * @param rootAtomElement l'élément racine XML dans lequel on souhaite lire un
	 * atome.
	 */
	xmlAtomParser(DOMElement * rootAtomElement);
	
	/**
	 * Destructeur de  la classe.
	 */
	virtual ~xmlAtomParser(void);
	
	/**
	 * Méthode pour construire un objet atom à partir des 
	 * informations présentes dans le document DOM.
	 *
	 * @return un objet de type atom qui a été initialisé à 
	 * partir des informations lues dans le document auquel 
	 * est attaché le parser.
	 */
	atom getAtom(void) const;	
};

#endif


