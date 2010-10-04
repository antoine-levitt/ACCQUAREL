/* 
 * The coniguration manager library of A.S.P.I.C. 
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
#ifndef _XML_ASPIC_CONFIGURATION_PARSER_
#define _XML_ASPIC_CONFIGURATION_PARSER_

#include <xmlDocumentParser.h>

/**
 * Class that contains all the informations about the basis sets
 * available.
 */
class basisDataBaseConfiguration
{
public:

	/**
	 * The relative path to the directory containing 
	 * all the basis XML files.
	 */
	string BasisDataBaseRelativePath;

	/**
	 * The relative path to the XML Schema for the
	 * basis files.
	 */
	string BasisDataBaseSchemaRelativePath;

	/**
	 * The relative path to the file containing
	 * the list of all the basis available.
	 */
	string BasisKeyListRelativePath;

	/**
	 * The relative path to the file containing
	 * the XML Schema for the list of all the basis available.
	 */
	string BasisKeyListSchemaRelativePath;

};

/**
 * Class that contains all the informations about the chemical elements
 * data base configuration
 */
class chemicalElementsDataBaseConfiguration
{
private:
	/**
	 * The relative path to the XML Document containing 
	 * the chemical elements data base.
	 */
	string ChemicalElementsDataBaseRelativePath;

	/**
	 * The relative path to the XML Schema for the
	 * chemical elements data base.
	 */
	string ChemicalElementsDataBaseSchemaRelativePath;

public:
	
	/**
	 * Default constructor.
	 */
	inline chemicalElementsDataBaseConfiguration(void)
		: ChemicalElementsDataBaseRelativePath("") ,
		ChemicalElementsDataBaseSchemaRelativePath("")
	{
		;
	}

	/**
	 * Default constructor.
	 */
	inline chemicalElementsDataBaseConfiguration(const string & chemicalElementsDataBaseRelativePath, const string & chemicalElementsDataBaseSchemaRelativePath)
		: ChemicalElementsDataBaseRelativePath(chemicalElementsDataBaseRelativePath) ,
		ChemicalElementsDataBaseSchemaRelativePath(chemicalElementsDataBaseSchemaRelativePath)
	{
		;
	}

	/**
	 * GET Method for the relative path to the chemical elements data base.
	 */
	inline const string & getChemicalElementsDataBaseRelativePath(void) const
	{
		return ChemicalElementsDataBaseRelativePath;
	}

	/**
	 * GET Method for the relative path to the chemical elements data base Schema.
	 */
	inline const string & getChemicalElementsDataBaseSchemaRelativePath(void) const
	{
		return ChemicalElementsDataBaseSchemaRelativePath;
	}

	/**
	 * SET Method for the relative path to the chemical elements data base.
	 */
	inline void setChemicalElementsDataBaseRelativePath(const string & chemicalElementsDataBaseRelativePath)
	{
		ChemicalElementsDataBaseRelativePath = chemicalElementsDataBaseRelativePath;
	}

	/**
	 * SET Method for the relative path to the chemical elements data base.
	 */
	inline void setChemicalElementsDataBaseSchemaRelativePath(const string & chemicalElementsDataBaseSchemaRelativePath)
	{
		ChemicalElementsDataBaseSchemaRelativePath = chemicalElementsDataBaseSchemaRelativePath;
	}

};

/**
 * Class that contains all the configurations informations about
 * the shell Types data base.
 *
 * This class contains the relative path to the XML document
 * that contains all the shell types.
 *
 * This class contains the relative path to the XML schema for
 * the shell types data base.
 */
class shellTypesDataBaseConfiguration
{
private:

	/**
	 * the string that contains the relative path to the 
	 * XML document that contains the shell types data base.
	 */
	string ShellTypesDataBaseRelativePath;

	/**
	 * The relative path to the XML Schema for the
	 * basis files.
	 */
	string ShellTypesDataBaseSchemaRelativePath;

public:

	/**
	 * Default constructor.
	 *
	 * This constructor build an empty object.
	 */
	inline shellTypesDataBaseConfiguration(void)
		: ShellTypesDataBaseRelativePath("") , ShellTypesDataBaseSchemaRelativePath("")
	{
		;
	}

	/**
	 * Constructor with specification.
	 *
	 * This constructor build an object which contains .
	 */
	inline shellTypesDataBaseConfiguration(const string & shellTypesDataBaseRelativePath , const string & shellTypesDataBaseSchemaRelativePath)
		: ShellTypesDataBaseRelativePath(shellTypesDataBaseRelativePath) , 
		ShellTypesDataBaseSchemaRelativePath(shellTypesDataBaseSchemaRelativePath)
	{
		;
	}

	/**
	 * Method GET for the relative path to the 
	 * XML document containing the shell type data base.
	 */
	inline const string & getShellTypesDataBaseRelativePath(void) const
	{
		return ShellTypesDataBaseSchemaRelativePath;
	}

	/**
	 * Method GET for the relative path to the 
	 * XML Schema for the shell types data base.
	 */
	inline const string & getShellTypesDataBaseSchemaRelativePath(void) const
	{
		return ShellTypesDataBaseSchemaRelativePath;
	}

	/**
	 * Method SET for the relative path to the 
	 * XML document containing the shell type data base.
	 */
	inline void setShellTypesDataBaseRelativePath(const string & shellTypesDataBaseSchemaRelativePath)
	{
		ShellTypesDataBaseSchemaRelativePath = shellTypesDataBaseSchemaRelativePath;
	}

	/**
	 * Method SET for the relative path to the 
	 * XML Schema for the shell types data base.
	 */
	inline void setShellTypesDataBaseSchemaRelativePath(const string & shellTypesDataBaseSchemaRelativePath)
	{
		ShellTypesDataBaseSchemaRelativePath = shellTypesDataBaseSchemaRelativePath;
	}

};


/**
 * class that parse the XML document $(ASPICROOT)/data/aspicConfiguration.xml.
 *
 * The purpose of this class is to read the user options for the A.S.P.I.C plate
 * forme. The options are the following :
 *
 * - Shall A.S.P.I.C use local datas ? (yes or no) : this is use to tell if
 * local file should be use for the gaussian basis, shell types, and chemical 
 * elements.
 *
 * - A.S.P.I.C server URI (string) : this string contains the base URI of
 * the server that hosts the data bases that will be used by the A.S.P.I.C plate 
 * form. 
 * This varaible contains something like file:///path/to/another/A.S.P.I.C/ (if the 
 * files are on the same computer) or something like 
 * http://cermics.enpc.fr/~lodier/A.S.P.I.C if the files are located on another computer.
 * Let us remark that this variable will not be used if the option use local datas is enabled.
 *
 * - The configuration of the chemical elements data base :
 * This falls into two parts, the URI of the XML document containing the data base
 * and the URI of the XML schema that validate the document.
 *
 * - The configuration of the shell types data base :
 * This falls into two parts, first the URI of the XML document containing the data
 * base, second the URI of the XML Schema that validates the XML document containing
 * the data base.
 *
 * - The configuration of the gaussian basis data base :
 *
 */
class xmlAspicConfigurationParser : public xmlDocumentParser
{
private:

	xmlAspicConfigurationParser * LocalAspicConfiguration;
	xmlAspicConfigurationParser * ExternalAspicConfiguration;

protected:

	/**
	 * Default and only constructor.
	 *
	 * This method constructs an object that will be abble to read the configuration 
	 * file of the A.S.P.I.C plate forme from the file 
	 * $(ASPICROOT)/data/aspicConfiguration.xml
	 */
	xmlAspicConfigurationParser(void);

	/**
	 * Destructor of class.
	 */
	virtual ~xmlAspicConfigurationParser(void);

	/**
	 * Method that tells if we should use the local data bases
	 * or a delocated data bases set.
	 *
	 * @return true if the local data bases are in use, false
	 * if not.
	 */
	bool _useLocalDataBases(void) const;

	/**
	 * Method that tells which data bases server shall be used.
	 *
	 * @retur the URI of the data base server for the A.S.P.I.C
	 * server.
	 *
	 * @warning If the configuration uses local data bases an empty string
	 * is returned.
	 */
	const string _getAspicServerRoot(void) const;

	/**
	 * Method that reads, in the configuration file, the configuration for the
	 * basis configuration.
	 *
	 * @return all the informations that can be found in the configuration
	 * file for the basis sets.
	 */
	const basisDataBaseConfiguration _getBasisDataBaseConfiguration(void) const;

	/**
	 * Method that reads, in the configuration file, the relative
	 * path to the chemical elements data base and the schema.
	 *
	 *
	 * 
	 * @return the relative path to the chemical elements data base.
	 */
	const chemicalElementsDataBaseConfiguration _getChemicalElementsDataBaseConfiguration(void) const;

	/**
	 * Method that reads, in the configuration file, the relative
	 * path to the shell types data base and the schema.
	 *
	 *
	 * 
	 * @return the relative path to the shell types data base and the XML schema.
	 */
	const shellTypesDataBaseConfiguration _getShellTypesDataBaseConfiguration(void) const;

public:

	/**
	 * This method is used to know if the plate forme should use local
	 * data bases. 
	 *
	 * - If the method returns true, it implies all data are gonna be read
	 * in the $(ASPICROOT) directory. 
	 * - If the method returns false, the data will be fetch from another 
	 * repository.
	 * 
	 * @return true if the option is to use local data bases, false if not.
	 */
	static bool useLocalDataBases(void);

	/**
	 * This method returns the base URI of the data bases repository.
	 *
	 * This value as sense if and only if the configuration uses delocated
	 * data bases. If this is the case an empty string will be returned.
	 *
	 * @return the URI of the ASPIC data base server.
	 */
	static const string getAspicServerRoot(void);

};

#endif