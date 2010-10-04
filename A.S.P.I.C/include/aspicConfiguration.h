/* 
 * The coniguration manager library of A.S.P.I.C. 
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
#ifndef _ASPIC_CONFIGURATION_
#define _ASPIC_CONFIGURATION_

#include <containor.h>
#include <iostream>
#include <string>
#include <xmlDocumentParser.h>
using namespace std;

class aspicConfiguration 
{
private:

protected:

public:

	/**
	 * This method is use to find the root 
	 * of the ASPIC files.
	 * 
	 * When the ASPIC_ROOT environment variable
	 * is set, then the value stored in this varible
	 * is returned.
	 *
	 * When the ASPIC_ROOT environment varible is
	 * not set then we try to look in different places :
	 * first "/etc/aspic"
	 * then  "/usr/local/aspic".
	 *
	 * If nothing was found this print an error message
	 * and the calling program will be aborted.
	 */
	static const string getAspicRoot(void);

	/**
	 * Méthode pour retrouver le chemin vers le fichier XML qui contient la base de données des 
	 * éléments chimiques.
	 *
	 * Lorsque le fichier XML ne peut pas être trouvé alors, une chaine 
	 * vide est renvoyée.
	 */
	static const string getChemicalElementsDataBasePath(void);
	
	/**
	 * Méthode pour retrouver le chemin vers le fichier schéma de la 
	 * base de données des éléments chimiques.
	 *
	 * Lorsque le fichier schéma ne peut pas être trouvé alors, une chaine 
	 * vide est renvoyée.
	 */
	static const string getChemicalElementsDataBaseSchemaPath(void);
	
	/**
	 * Méthode pour retrouver le répertoire racine qui contient les
	 * tout les fichiers de bases gaussiennes.
	 */
	static const string getBasisDataBaseRoot(void);

	/**
	 * Méthode pour retrouver le chemin vers le fichier XML de la 
	 * base de données des bases chimiques.
	 *
	 * Lorsque le fichier XML ne peut pas être trouvé alors, une chaine vide est renvoyée.
	 */
	static const string getBasisDataBasePath(const string & basisName);

	/**
	 * Méthode pour retrouver le schéma des fichiers XMLs qui contiennent les
	 * base gaussiennes.
	 *
	 * Lorsque le fichier schéma n'a pas été trouvé alors on renvoie une chaine vide.
	 */
	static const string getBasisDataBaseSchemaPath(void);
	
	/**
	 * Méthode pour retrouver le schéma des fichiers XMLs qui contiennent une fonction de base.
	 *
	 * Lorsque le fichier schéma n'a pas été trouvé alors on renvoie une chaine vide.
	 */
	static const string getBasisFunctionSchemaPath(void);
	
	/**
	 * Méthode pour retrouver la liste de toute les bases gaussiennes existantes.
	 */
	static const string getBasisKeyListPath(void);
		
	/**
	 * Méthode pour retrouver la liste de toute les bases gaussiennes existantes.
	 */
	static const string getBasisKeyListSchemaPath(void);

	
	/**
	 * Méthode pour retrouver le fichier schéma pour les molécules.
	 *
	 * Lorsque le fichier XML n'est pas trouvé une chaine vide est renvoyée.
	 */
	static const string getMoleculeSchemaPath(void);

	/**
	 * Méthode pour retrouver le chemin vers le fichier XML de la 
	 * base de données des éléments chimiques.
	 *
	 * Lorsque le fichier XML ne peut pas être trouvé alors, une chaine vide est renvoyée.
	 */
	static const string getShellTypesDataBasePath(void);

	/**
	 * Méthode pour retrouver le chemin vers le fichier schéma de la 
	 * base de données des éléments chimiques.
	 *
	 * Lorsque le fichier schéma ne peut pas être trouvé alors, une chaine 
	 * vide est renvoyée.
	 */
	static const string getShellTypesDataBaseSchemaPath(void);	

	/**
	 * Méthode pour retrouver le chemin vers le fichier schéma d'une matrice
	 *
	 * Lorsque le fichier schéma ne peut pas être trouvé alors, une chaine 
	 * vide est renvoyée.
	 */
	static const string getMatrixSchemaPath(void);	
};


/**
 * Function to know if a directory exists.
 * 
 * This function returns true if the string "path" 
 * given as an argument is a directory. If the is a file
 * or does not exists it will return false.
 */
extern bool isDirectory(const string & path);

/**
 * Fonction pour savoir si un fichier existe et est accessible en lecture.
 *
 * @param path le chemin vers le fichier.
 *
 * @return vrai lorsque le fichier est accessible en lecture, faux lorsque le fichier n'existe pas
 * ou si l'utilsateur du programme ne possède pas les droit pour lire le fichier.
 */
extern bool isFileReadable(const string & path);

#endif

