/* 
* The containor library of A.S.P.I.C. 
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
#ifndef _XML_IPOINT_3_PARSER_
#define _XML_IPOINT_3_PARSER_

#include "ipoint.h"
#include "xmlParser.h"


/**
 * Classe pour lire un triplet d'entier dans un document XML.
 */
class xmlIPoint3Parser : public virtual xmlParser
{

private:
	string IPoint3TagName;
	string XTagName;
	string YTagName;
	string ZTagName;

public :
	/**
	 * Constructeur.
	 */
	xmlIPoint3Parser(DOMElement * DOMPointElement);

	/**
	 * Constructeur.
	 */
	xmlIPoint3Parser(void);

	/**
	 * Destructeur.
	 */
	virtual ~xmlIPoint3Parser(void);

	/**
	 * Méthode pour accéder à la valeur d'un triplet d'entier
	 * dans un element DOM qui a pour racine rootElement.
	 */
	ipoint<3> getIPoint3(const DOMElement * rootElement) const;

	/**
	 * Méthode pour accéder à la valeur d'un triplet d'entier
	 * contenu dans le document.
	 */
	ipoint<3> getIPoint3(void) const;

	/**
	 * Méthode pour accéder au tag qui contient le triplet d'entier.
	 */
	const string & getIPoint3TagName(void) const;

	/**
	 * Méthode pour accéder au tag qui contient la première coordonées du triplet d'entiers.
	 */
	const string & getXTagName(void) const;

	/**
	 * Méthode pour accéder au tag qui contient la seconde coordonées du triplet d'entiers.
	 */
	const string & getYTagName(void) const;

	/**
	 * Méthode pour accéder au tag qui contient la triosième coordonées du triplet d'entiers.
	 */
	const string & getZTagName(void) const;
	
	/**
	 * Méthode pour accéder au tag qui contient le triplet d'entier.
	 */
	void setIPoint3TagName(const string & iPoint3DTagName);

	/**
	 * Méthode pour accéder au tag qui contient la première coordonées du triplet d'entiers.
	 */
	void setXTagName(const string & xTagName);

	/**
	 * Méthode pour accéder au tag qui contient la seconde coordonées du triplet d'entiers.
	 */
	void setYTagName(const string & yTagName);

	/**
	 * Méthode pour accéder au tag qui contient la triosième coordonées du triplet d'entiers.
	 */
	void setZTagName(const string & zTagName);

	/**
	 * Méthode pour remettre la valeur des tag à leur valeur par défaut.
	 */
	void resetTagNames(void);
};


#endif

