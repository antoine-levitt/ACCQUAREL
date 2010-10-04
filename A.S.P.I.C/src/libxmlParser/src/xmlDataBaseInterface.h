/* 
 * The xml parser library of A.S.P.I.C. 
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
#ifndef _XML_DATA_BASE_INTERFACE_
#define _XML_DATA_BASE_INTERFACE_

#include "xmlDocumentParser.h"

/**
 * Generic class to manage XML data base.
 */
class xmlDataBaseInterface : protected xmlDocumentParser
{
private:

	/**
	 * Le nom du tag qui contient une entrée.
	 */
	const string EntryTagName;
	
	/**
	 * Le nom de l'attribut qui est la clé de l'élément.
	 */
	const string IdAttributeName;
	
protected:
	
	/**
	 * Chaine de charactère qui contient le chemin  vers le fichier de la base de données.
	 */
	virtual string getDocumentURI(void) const =0;
	
	/**
	 * Méthode pour récupérer uniquement la partie
	 * du document qui contient les informations
	 * sur l'élément que nous recherchons.
	 */
	DOMElement * getElementById(const string & key) const;
	
	/**
	 * Méthode GET pour accéder à une entrée de la base de données.
	 */
	DOMNode * getDataBaseEntry(int item) const;
	
	/**
	 * Method GET for the ID of an entry in the data base.
	 */
	string getID4Entry(int item) const;
	
	/**
	 * Chaine de charactére qui contient le chemin vers le fichier de schéma.
	 */
	virtual string getSchemaURI(void) const =0;

	/**
	 * Méthode GET pour le nombre d'entrées dans la base de données
	 */
	int getNbrOfEntries(void) const;
	
public:
	
	/**
	 * Le constructeur.
	 */
		xmlDataBaseInterface(const string & entryTagName , const string & idAttributeName);
	
	/**
		* Le destructeur.
	 */
	virtual ~xmlDataBaseInterface(void);
};



#endif

