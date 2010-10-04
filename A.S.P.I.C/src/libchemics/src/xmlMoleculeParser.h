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
#ifndef _XML_MOLECULE_PARSER_
#define _XML_MOLECULE_PARSER_

#include "atom.h"
#include "xmlBasisDataBaseInterface.h"
#include "molecule.h"
#include "xmlAtomParser.h"
#include <xmlParser.h>


class xmlMoleculeParser : virtual public xmlParser
{
private:

protected:
	/**
	 * Method GET to access a node of the doccument 
	 * that correspnds to the atom number item.
	 */
	DOMElement * getAtomDOMElement(int item) const;

public:

	/**
	 * Constructor.
	 */
	xmlMoleculeParser(void);

	/**
	 * Constructor.
	 */
	xmlMoleculeParser(DOMElement * rootElement);

	/**
	 * Destructor.
	 */
	virtual ~xmlMoleculeParser(void);
	
	/**
	 * String that contains the tag name for atoms in the molecule.
	 */
	static const string getAtomTagName(void);
	
	/**
	 * String that contains the tag name for atom basis.
	 */
	static const string getBasisNameTagName(void);
	
	/**
	 * String that contains the tag name for the charge of the molecule.
	 */
	static const string getChargeTagName(void);

	/**
	 * Méthode pour connaitre le nom du tag du document qui contient le 
	 * description de la molécule.
	 *
	 * Dans le document xml une description de la molécule est disponible,
	 * cette description est dans un élément dont le nom de tag est renvoyé par 
	 * cette méthode.
	 *
	 * @return le nom du tag qui contient la description de la molécule.
	 */
	static const string getDescriptionTagName(void);
		
	/**
	 * Method GET for the whole molecule.
	 */
	molecule getMolecule(void) const;

	/**
   * String that contains the tag name for the molecule name.
	 */
	static const string getNameTagName(void);

	/**
	 * Méthode pour connaitre le nom des balises qui encadre la
	 * liste des éléments (atome + base) qui composent la molécule
	 */
	static const string getElementsListTagName(void);

	/**
	 * Méthode pour connaitre le nom des balises qui contiennent la description
	 * d'un élément (atome + base) qui compose la molécule.
	 */
	static const string getElementTagName(void);

	/**
	 * Méthode pour connaitre le nom des balises qui lcontiennent la clé de 
	 * de la base de l'élément qui compose la molécule.
	 */
	static const string getElementBasisKeyTagName(void);

};

#endif

