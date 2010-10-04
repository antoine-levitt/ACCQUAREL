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
#ifndef _XML_GAUSSIAN_BASIS_SHELL_TYPE_PARSER_
#define _XML_GAUSSIAN_BASIS_SHELL_TYPE_PARSER_

#include "shellType.h"
#include <xmlParser.h>


/**
 * Classe pour parser un DOMElement (morceau de fichier XML) qui contient un
 * type de couche.
 * Le propos de cette classe est de construire un objet shellType contenant
 * les inforamtions qui sont présente dans le document.
 */
class xmlShellTypeParser : public virtual xmlParser
{
private:

protected:

public:
	
	/**
	 * Le constructeur par défaut.
	 */
	xmlShellTypeParser(void);

	/**
	 * Le constructeur avec la racine du document à parser.
	 */
	xmlShellTypeParser(DOMElement * rootShellTypeElement);

	/**
	 * Le destructeur.
	 */
	virtual ~xmlShellTypeParser(void);

	/**
   * Méthode qui lit et construit l'objet shellType présent dans le document.
	 */
	shellType getShellType(void) const;

	/**
   * Méthode qui lit et construit l'objet shellType présent dans le document.
	 */
	shellType getShellType(const DOMElement * rootElement) const;
	
	/**
	 * Le nom du tag qui contient le nom de la couche.
	 */
	static const string getShellTypeNameTagName(void);

	/**
	 * Le nom du tag qui contient le nom de la couche.
	 */
	static const string getShellTypeKeyAttributeName(void);

	/**
	 * Le nom du tag qui contient le le type de couche.
	 */
	static const string getShellTypeTagName(void);

	/**
	 * Le nom du tag qui contient le degré d'un monome présent dans la couche.
	 */
	static const string getMonomeDegreeTagName(void);

	/**
	 * Le nom du tag qui contient le degré d'un monome présent dans la couche.
	 */
	static const string getMonomeDegreesListTagName(void);

};

#endif

