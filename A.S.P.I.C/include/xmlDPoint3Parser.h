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
#ifndef _XML_DPOINT_3_PARSER_
#define _XML_DPOINT_3_PARSER_

#include "dpoint.h"
#include <xmlParser.h>

/**
 * Classe pour lire une position dans un document XML.
 *
 * Cette classe permet de lire des document de la forme :
 *
 * <Position>
 *		<x> 1 </x>
 *		<y> 2 </y> 
 *		<z> 3 </z>
 * </Position>
 *
 * A partir d'un tel document (ou partie de document) cette classe 
 * construit un objet dpoint<3> qui a pour valeur (1,2,3) si on l'applique
 * à l'exemple ici proposé.
 */
class xmlDPoint3Parser : virtual public xmlParser
{
private:
	
	string DPoint3TagName;
	string XTagName;
	string YTagName;
	string ZTagName;


public :
	
	/**
	 * Constructeur.
	 */
	xmlDPoint3Parser(void);

	/**
	 * Constructeur.
	 */
	xmlDPoint3Parser(DOMElement * DOMPointElement);

	/**
	 * Destructeur.
	 */
	virtual ~xmlDPoint3Parser(void);

	
	/**
	 * Méthode pour accéder à la valeur.
	 *
	 * @return la valeur de l'objet dpoint<3> lu dans le 
	 * document auquel le parser est attaché.
	 */
	dpoint<3> getDPoint3(const DOMElement * root) const;


	/**
	 * Méthode pour accéder à la valeur.
	 *
	 * @return la valeur de l'objet dpoint<3> lu dans le 
	 * document auquel le parser est attaché.
	 */
	dpoint<3> getDPoint3(void) const;

	/**
	 * Méthode pour accéder au nom de la balise qui contient 
	 * la positon.
	 */
	const string & getDPoint3TagName(void) const;

		/**
	 * Méthode pour accéder au nom de la balise qui contient 
	 * la première coordonnée de positon.
	 */
	const string & getXTagName(void) const;

		/**
	 * Méthode pour accéder au nom de la balise qui contient 
	 * la seconde coordonnées de positon.
	 */
	const string & getYTagName(void) const;

		/**
	 * Méthode pour accéder au nom de la balise qui contient 
	 * la troisième coordonnée positon.
	 */
	const string & getZTagName(void) const;

	/**
	 * Méthode pour accéder au nom de la balise qui contient 
	 * la positon.
	 */
	void setDPoint3TagName(const string & dPoint3TagName);

	/**
	 * Méthode pour accéder au nom de la balise qui contient 
	 * la première coordonnée de positon.
	 */
	void setXTagName(const string & xTagName);

	/**
	 * Méthode pour accéder au nom de la balise qui contient 
	 * la seconde coordonnées de positon.
	 */
	void setYTagName(const string & yTagName);

	/**
	 * Méthode pour accéder au nom de la balise qui contient 
	 * la troisième coordonnée positon.
	 */
	void setZTagName(const string & zTagName);

	/** 
	 * Méthode pour remettre les tags à leurs valeurs initiale.
	 */
	void resetTagNames(void);
};

#endif


