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
#include "aspicConfiguration.h"
#include <sys/stat.h>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to know if a directory exists.
// 
// This function returns true if the string "path" 
// given as an argument is a directory. If the is a file
// or does not exists it will return false.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool isDirectory(const string & path)
{
	
	struct stat fileStatus;

	stat(path.c_str() , &fileStatus);

	if(S_IFDIR == fileStatus.st_mode) 
		return true;
	
	else
		return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to know if a directory exists.
// 
// This function returns true if the string "path" 
// given as an argument is a directory. If the is a file
// or does not exists it will return false.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool isFileReadable(const string & path)
{
	if(path.empty()) {
		return false;
	}
	
	FILE * fileDescriptor = fopen(path.c_str() , "r");
	
	if(fileDescriptor == NULL) {
		return false;
	} else {
		fclose(fileDescriptor);
		return true;
	}	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// static const string aspicConfiguration::getAspicRoot(void) 
//
// This method is use to find the root 
// of the ASPIC files.
// 
// When the ASPICROOT environment variable
// is set, then the value stored in this varible
// is returned.
//
// When the ASPIC_ROOT environment varible is
// not set then we try to look in different places :
// first "/etc/aspic"
// then  "/usr/local/aspic".
//
// If nothing was found this print an error message
// and the calling program will be aborted.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getAspicRoot(void) 
{
	string aspicRoot;
	char * aspicRootCh = getenv("ASPICROOT");
	///////////////////////////////////////////////////////////////////////////////////////
	// 1 - The environement varible ASPIC_ROOT is set.
	if (aspicRootCh != NULL) {
		aspicRoot = aspicRootCh;
	} else {
	aspicRoot = "";
}

	// If the environment varible has a value we stop here ...
	if(! aspicRoot.empty()) {
		return aspicRoot;
	}	

	/////////////////////////////////////////////////////////
	// 2 - The directory "/etc/aspic".
	aspicRoot = "/etc/aspic";

	if(isDirectory(aspicRoot)) {
		return aspicRoot;
	}

	/////////////////////////////////////////////////////////
	// 3 - The directory "/usr/local/aspic".
	aspicRoot = "/usr/local/aspic";

	if(isDirectory(aspicRoot)) {
		return aspicRoot;
	}

	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getAspicRoot(void)" << endl;
	cerr << "Error : the program was unable to find the datas for the aspic." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure the directory \"/etc/aspic\" or \"/usr/local/aspic\" exists." << endl;
	cerr << "Error :" << endl;
	cerr << "Error : If you install aspic in a specific location make sure the " << endl;
	cerr << "Error : ASPICROOT environment varible is set." << endl;

	exit(1);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// const string & aspicConfiguration::getBasisDataBasePath(const string & basisName)
//
// Méthode pour retrouver le fichier XML associé à une base gaussienne.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getBasisDataBasePath(const string & basisName)
{
	string gaussianBasisDB = getBasisDataBaseRoot();
	gaussianBasisDB += "/";
	gaussianBasisDB += basisName;
	gaussianBasisDB +=".xml";
	
	if(isFileReadable(gaussianBasisDB)) {
		return gaussianBasisDB;
	}
	
	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getBasisDataBasePath(const string basisName)" << endl;
	cerr << "Error : unable to find the XML file for the gaussian basis DB." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure you correctly set the name of the basis :\"" << basisName << "\"" << endl;
	cerr << "Error : Make sure the file with path \"" << gaussianBasisDB << "\" exists." << endl;
	cerr << "Error : Make sure the ASPICROOT environement variable is corectly set." << endl;
	return "";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// const string & aspicConfiguration::getBasisDataBaseSchemaPath(void)
//
// Méthode pour retrouver le schéma des fonctions de base.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getBasisDataBaseSchemaPath(void)
{
	string gaussianBasisDBSchema = getAspicRoot();
	gaussianBasisDBSchema += "/data/chemics/gaussianBasis.xsd";
	
	if(isFileReadable(gaussianBasisDBSchema)) {
		return gaussianBasisDBSchema;
	}
	
	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getBasisDataBaseSchemaPath(void)" << endl;
	cerr << "Error : unable to find the XML file for the gaussian basis DB." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure the file with path \"" << gaussianBasisDBSchema << "\" exists." << endl;
	cerr << "Error : Make sure the ASPICROOT environement variable is corectly set." << endl;
	
	return "";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// const string & aspicConfiguration::getBasisFunctionSchemaPath(void)
//
// Méthode pour retrouver la liste de toutes les base présentes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getBasisFunctionSchemaPath(void)
{	
	string basisFunctionSchema = getAspicRoot();
	basisFunctionSchema += "/data/chemics/basisFunction.xsd";

	if(isFileReadable(basisFunctionSchema)) {
		return basisFunctionSchema;
	}

	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getBasisFunctionSchemaPath(void)" << endl;
	cerr << "Error : unable to find the XML file for the basis function." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure the file with path \"" << basisFunctionSchema << "\" exists." << endl;
	cerr << "Error : Make sure the ASPICROOT environement variable is corectly set." << endl;
	
	return "";
}
	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// const string & aspicConfiguration::getBasisListKeyPath(void)
//
// Méthode pour retrouver la liste de toutes les base présentes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getBasisKeyListPath(void)
{	
	string basisListFileName = getAspicRoot();
	basisListFileName += "/data/chemics/gaussianBasis/basisKeyList.xml";

	if(isFileReadable(basisListFileName)) {
		return basisListFileName;
	}
	
	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getBasisKeyListPath(void)" << endl;
	cerr << "Error : unable to find the XML file for the gaussian basis key list." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure the file with path \"" << basisListFileName << "\" exists." << endl;
	cerr << "Error : Make sure the ASPICROOT environement variable is corectly set." << endl;
	
	return "";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// const string & aspicConfiguration::getBasisListKeySchemaPath(void)
//
// Méthode pour retrouver la liste de toutes les base présentes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getBasisKeyListSchemaPath(void)
{	
	string basisListFileName = getAspicRoot();
	basisListFileName += "/data/chemics/basisKeyList.xsd";

	if(isFileReadable(basisListFileName)) {
		return basisListFileName;
	}
	
	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getBasisKeyListPath(void)" << endl;
	cerr << "Error : unable to find the schema of XML file for the gaussian basis key list." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure the file with path \"" << basisListFileName << "\" exists." << endl;
	cerr << "Error : Make sure the ASPICROOT environement variable is corectly set." << endl;
	
	return "";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// const string aspicConfiguration::getBasisDataBaseRoot(void)
//
// Méthode pour connaitre le répertoire qui contient tout les fichiers de base.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getBasisDataBaseRoot(void)
{
	string basisRoot = getAspicRoot();
	basisRoot += "/data/chemics/gaussianBasis";
	return basisRoot;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// static const string aspicConfiguration::getChemicalElementsDataBasePath(void) 
//
// Cette méthode permet de retourver le fichier XML qui contient la base de données des éléments chimiques.
// Normalement ce fichier devrait se trouver dans $(ASPICROOT)/data/chemicalElements/chemicalElements.xml.
// 
// On va faire les choses comme il faut, c'est à dire que nous allons tester ce chemin (existance du fichier)
// et s'il existe on renvoi le chemin trouvé, sinon on renvoie la chaine vide.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getChemicalElementsDataBasePath(void) 
{
	string chemicalElmentsDB = getAspicRoot();
	chemicalElmentsDB += "/data/chemics/chemicalElements.xml";
	
	if(isFileReadable(chemicalElmentsDB)) {
		return chemicalElmentsDB;
	}
	
	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getChemicalElementsDataBasePath(void)" << endl;
	cerr << "Error : unable to find the XML file for the chemical elements DB." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure the file with path \"" << chemicalElmentsDB << "\" exists." << endl;
	cerr << "Error : Make sure the ASPICROOT environement variable is corectly set." << endl;
	
	return "";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// static const string aspicConfiguration::getChemicalElementsDataBaseSchema(void) 
//
// Cette méthode permet de retourver le fichier schéma qui décrit la base de données des éléments chimiques.
// Normalement ce fichier devrait se trouver dans $(ASPICROOT)/data/chemics/chemicalElements.xsd.
// 
// On va faire les choses comme il faut, c'est à dire que nous allons tester ce chemin (existance du fichier)
// et s'il existe on renvoi le chemin trouvé, sinon on renvoie la chaine vide.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getChemicalElementsDataBaseSchemaPath(void) 
{
	string chemicalElmentsDBSchema = getAspicRoot();
	chemicalElmentsDBSchema += "/data/chemics/chemicalElements.xsd";
	
	if(isFileReadable(chemicalElmentsDBSchema)) {
		return chemicalElmentsDBSchema;
	}
	
	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getChemicalElementsDataBaseSchema(void)" << endl;
	cerr << "Error : unable to find the schema file for the chemical elements DB." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure the file with path \"" << chemicalElmentsDBSchema << "\" exists." << endl;
	cerr << "Error : Make sure the ASPICROOT environement variable is corectly set." << endl;
	
	return "";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// const string aspicConfiguration::getMoleculeSchemaPath(void)
//
// Méthode pour le fichier de schéma pour les molécules.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getMoleculeSchemaPath(void)
{	
	string schema = getAspicRoot();
	schema+= "/data/chemics/molecule.xsd";

	if(isFileReadable(schema)) {
		return schema;
	}
	
	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getMoleculeSchemaPath(void)" << endl;
	cerr << "Error : unable to find the schema of XML file for the molecule." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure the file with path \"" << schema << "\" exists." << endl;
	cerr << "Error : Make sure the ASPICROOT environement variable is corectly set." << endl;
	
	return "";
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// static const string aspicConfiguration::getChemicalElementsDataBasePath(void) 
//
// Cette méthode permet de retourver le fichier XML qui contient la base de données des éléments chimiques.
// Normalement ce fichier devrait se trouver dans $(ASPICROOT)/data/chemics/shellTypes.xml.
// 
// On va faire les choses comme il faut, c'est à dire que nous allons tester ce chemin (existance du fichier)
// et s'il existe on renvoi le chemin trouvé, sinon on renvoie la chaine vide.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getShellTypesDataBasePath(void) 
{
	string shellTypesDB = getAspicRoot();
	shellTypesDB += "/data/chemics/shellTypes.xml";
	
	if(isFileReadable(shellTypesDB)) {
		return shellTypesDB;
	}
	
	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getShellTypesDataBasePath(void)" << endl;
	cerr << "Error : unable to find the XML file for the shell types DB." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure the file with path \"" << shellTypesDB << "\" exists." << endl;
	cerr << "Error : Make sure the ASPICROOT environement variable is corectly set." << endl;
	
	return "";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// static const string aspicConfiguration::getChemicalElementsDataBaseSchema(void) 
//
// Cette méthode permet de retourver le fichier schéma qui décrit la base de données des éléments chimiques.
// Normalement ce fichier devrait se trouver dans $(ASPICROOT)/data/chemics/shellTypes.xsd.
// 
// On va faire les choses comme il faut, c'est à dire que nous allons tester ce chemin (existance du fichier)
// et s'il existe on renvoi le chemin trouvé, sinon on renvoie la chaine vide.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getShellTypesDataBaseSchemaPath(void) 
{
	string shellTypesDBSchema = getAspicRoot();
	shellTypesDBSchema += "/data/chemics/shellTypes.xsd";
	
	if(isFileReadable(shellTypesDBSchema)) {
		return shellTypesDBSchema;
	}
	
	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getShellTypesDataBaseSchema(void)" << endl;
	cerr << "Error : unable to find the schema file for the shell types DB." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure the file with path \"" << shellTypesDBSchema << "\" exists." << endl;
	cerr << "Error : Make sure the ASPICROOT environement variable is corectly set." << endl;
	
	return "";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// static const string aspicConfiguration::getChemicalElementsDataBaseSchema(void) 
//
// Cette méthode permet de retourver le fichier schéma qui décrit les matrices.
// Normalement ce fichier devrait se trouver dans $(ASPICROOT)/data/matrix/matrix.xsd.
// 
// On va faire les choses comme il faut, c'est à dire que nous allons tester ce chemin (existance du fichier)
// et s'il existe on renvoi le chemin trouvé, sinon on renvoie la chaine vide.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string aspicConfiguration::getMatrixSchemaPath(void) 
{
	string matrixSchema = getAspicRoot();
	matrixSchema += "/data/matrix/matrix.xsd";
	
	if(isFileReadable(matrixSchema)) {
		return matrixSchema;
	}
	
	////////////////////////////////////////////////////////////
	// Error : unable to find aspic root. We write a message and
	// abort the program.
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error : const string in aspicConfiguration::getMatrixSchemaPath(void)" << endl;
	cerr << "Error : unable to find the schema file for the shell types DB." << endl;
	cerr << "Error : ---------------------------------------------------------------" << endl;
	cerr << "Error :" << endl;
	cerr << "Error : Make sure the file with path \"" << matrixSchema << "\" exists." << endl;
	cerr << "Error : Make sure the ASPICROOT environement variable is corectly set." << endl;
	
	return "";
}


