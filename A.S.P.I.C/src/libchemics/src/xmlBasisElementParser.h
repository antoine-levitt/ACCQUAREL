/* 
* The chemics library of A.S.P.I.C. 
 * Written and directed by Fran�ois Lodier support.aspic@gmail.com.
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

#include <xmlParser.h>
#include "basisElement.h"
/**
 * Classe dont le propos est de transofrmer une partie de document
 * XML en objet basisElement.
 */
class xmlBasisElementParser : public xmlParser
{
private:
protected:
public:
	/**
	 * Constructeur par d�faut.
	 */
	xmlBasisElementParser(void);
	
	/**
	 * Constructeur avec un �l�ment racine.
	 */
	xmlBasisElementParser(DOMElement * rootElement);
	
	/**
	 * Destructeur.
	 */
	virtual ~xmlBasisElementParser(void);
	
	/**
	 * Méthode pour récupérer l'élément de base.
	 */
	basisElement getBasisElement(void) const;

	/**
	 * Méthode pour r�cup�rer l'élément de base.
	 */
	basisElement getBasisElement(const DOMElement * rootElement) const;

	/**
	 * Méthode pour récupérer le nom du tag qui contient un élément de base.
	 */
	static const string getBasisElementTagName(void);
	
	/**
	 * Méthode pour récupérer le nom du tag qui contient un élément de base.
	 */
	static const string getElementKeyAttributeName(void);

	/**
	 * Méthode pour récupérer le nom du tag qui contient un élément de base.
	 */
	static const string getBasisKeyAttributeName(void);

	/**
	 * Méthode pour récupérer le nom du tag qui contient une couche d'élément de base.
	 */
	static const string getShellTagName(void);
	
	/**
	 * Méthode pour récupérer le nom du tag qui contient une couche d'élément de base.
	 */
	static const string getShellsListTagName(void);
	
	/**
	 * Méthode pour récupérer le nom du tag qui contient le type de couche d'une couche d'élément de base.
	 */
	static const string getShellTypeKeyTagName(void);
	
	
};
