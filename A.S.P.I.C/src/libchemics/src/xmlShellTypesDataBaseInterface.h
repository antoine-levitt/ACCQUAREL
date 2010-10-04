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
#ifndef _GAUSSIAN_BASIS_SHELL_TYPE_DATA_BASE_INTERFACE_
#define _GAUSSIAN_BASIS_SHELL_TYPE_DATA_BASE_INTERFACE_

#include <xmlDataBaseInterface.h>
#include "shellType.h"
#include "xmlShellTypeParser.h"

class xmlShellTypesDataBaseInterface : public xmlDataBaseInterface
{
private :
	
	/**
	 * Le nom du tag qui contient une entrée dans la base de données.
	 */
	static const string	ShellTypeTagName;

	/**
	 * Le nom de l'attribut qui contient la clé de chaque élément dans la base de donnée.
	 */
	static const string ShellTypeKeyAttributeName;
	
	/**
	 * L'objet xmlShellTypeDataBase instancié lors d'une connection.
	 */
	static xmlShellTypesDataBaseInterface * ShellTypesDataBase;


	protected:
	
	/**
	 * Method GET to find the document location.
	 */
	virtual string getDocumentURI(void) const;
	
	/**
	 * Methog GET to find the schema location.
	 */
	virtual string getSchemaURI(void) const;

	/**
 	 * The constructor.
	 */
	xmlShellTypesDataBaseInterface(void);
	
	/**
 	 * The destructor.
	 */
	virtual ~xmlShellTypesDataBaseInterface(void);

	public:
	
	/**
	 * Méthode pour se connecter à la base de données.
	 */		
	static void connect(void);
	
	/**
	 * Méthode pour savoir si l'on est connecté à la base de données.
	 */		
	static bool connected(void);
	
	/**
	 * Méthode pour se deconnecter de la base de données.
	 */
	static void disconnect(void);
	
	/**
	 * Méthode pour connaitre le nom de l'attribut qui contient la clé pour les élément de la base de donnée.
	 */
	static const string & getShellTypeKeyAttibuteName(void);

	/**
	 * Méthode pour connaitre le nom de l'attribut qui contient la clé pour les élément de la base de donnée.
	 */
	static const string & getShellTypeTagName(void);

	/**
	 * Méthode pour connaitre le nombre de type de couches stockés dans la base de données 
	 */
	static int getNbrOfShellTypes(void);

	/**
	 * Méthode pour retrouver un type de couche à partir de sa clé.
	 */
	static shellType getShellType(const string & shellTypeKey);

	/**
	 * Méthode pour accéder à un élémént item de la base de donnée.
	 *
	 * @warning il faut que item soit positif ou nul et strictement inférieur
	 * aux nombre d'éléments contenus dans la base de donnée.
	 */
	static shellType getShellType(const int & item);
};


#endif


