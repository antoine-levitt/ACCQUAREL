/* 
 * The configuration management library of A.S.P.I.C. 
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
#include "xmlAspicConfigurationParser.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// The default and only constructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlAspicConfigurationParser::xmlAspicConfigurationParser(void)
: xmlDocumentParser()
{
	;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The destructor.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
xmlAspicConfigurationParser::~xmlAspicConfigurationParser(void)
{
	;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to know if local configuration should be used.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool xmlAspicConfigurationParser::_useLocalDataBases(void) const
{
	const string useLocalDataBaseTagName = "DataBasesType";
	DOMNode * useLocalDataBasesElement = getElementByTagName(useLocalDataBaseTagName);
	string useLocalDataBasesValue;

	if(useLocalDataBasesElement == NULL) {
		cerr << "Error : in bool xmlAspicConfigurationParser::_useLocalDataBases(void) const" << endl;
		cerr << "Error : no tag with name \"" << useLocalDataBaseTagName << "\" was found." << endl;
		cerr << "Error : unable to find data bases type for A.S.P.I.C." << endl;
		exit(1);
	}

	useLocalDataBasesValue = getNodeStringValue(useLocalDataBasesElement , xmlParser::Remove_White_Space);
	
	if(useLocalDataBasesValue == "Ext") {
		return false;
	}

	if(useLocalDataBasesValue == "Local") {
		return false;
	}
	
	cerr << "Error : in bool xmlAspicConfigurationParser::_useLocalDataBases(void) const" << endl;
	cerr << "Error : unknown data base type \"" << useLocalDataBasesValue << "\"." << endl;
	cerr << "Error : unable to find data bases type for A.S.P.I.C." << endl;
	exit(1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to find the base URI of the server for ASPIC databases.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string xmlAspicConfigurationParser::_getAspicServerRoot(void) const
{
	
	// If local data bases are in use there is no server.
	// so we just return an empty string.
	if(_useLocalDataBases()) {
		return "";
	}

	const string dataBaseServerURITagName = "DataBaseServerURI";
	DOMNode * dataBasesServerURIElement = getElementByTagName(dataBaseServerURITagName);
	string dataBasesServerURIValue;

	if(dataBasesServerURIElement == NULL) {
		cerr << "Error : in const string xmlAspicConfigurationParser::_getAspicServerRoot(void) const" << endl;
		cerr << "Error : no tag with name \"" << dataBaseServerURITagName << "\" was found." << endl;
		cerr << "Error : unable to find the data bases server for A.S.P.I.C." << endl;
		exit(1);
	}

	dataBasesServerURIValue = getNodeStringValue(dataBasesServerURIElement , xmlParser::Remove_White_Space);	

	if(dataBasesServerURIValue.empty() || dataBasesServerURIValue == "") {
		cerr << "Error : in const string xmlAspicConfigurationParser::_getAspicServerRoot(void) const" << endl;
		cerr << "Error : empty URI for data bases server was found." << endl;
		cerr << "Error : unable to find the data bases server for A.S.P.I.C." << endl;
		exit(1);	
	}
	return dataBasesServerURIValue;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to find the basis data base configuration.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
const basisDataBaseConfiguration xmlAspicConfigurationParser::_getBasisDataBaseConfiguration(void) const
{

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to find the chemical elements data base configuration.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
const chemicalElementsDataBaseConfiguration xmlAspicConfigurationParser::_getChemicalElementsDataBaseConfiguration (void) const
{	
	const string chemicalElementsDataBaseTagName = "ChemicalElementsDataBase";
	const string chemicalElementsDataBasePathTagName = "ChemicalElementsDataBasePath";
	const string chemicalElementsDataBaseSchemaPathTagName = "ChemicalElementsDataBaseSchemaPath";

	chemicalElementsDataBaseConfiguration chemicalElementsConfiguration("","");

	DOMElement * chemicalElementsDataBaseElement = (DOMElement *)getElementByTagName(chemicalElementsDataBaseTagName);
	
	if( chemicalElementsDataBaseElement== NULL) {
		cerr << "Error : in const pair<string,string> xmlAspicConfigurationParser::_getChemicalElementsDataBaseRelativePath(void) const" << endl;
		cerr << "Error : no tag with name \"" << chemicalElementsDataBaseTagName << "\" was found." << endl;
		cerr << "Error : no configuration for the chemical elements data base was found" << endl;
		return chemicalElementsConfiguration;
	}

	DOMNode * chemicalElementsDataBasePathElement = getElementByTagName(chemicalElementsDataBaseElement,chemicalElementsDataBasePathTagName);

	if( chemicalElementsDataBasePathElement== NULL) {
		cerr << "Error : in const pair<string,string> xmlAspicConfigurationParser::_getChemicalElementsDataBaseRelativePath(void) const" << endl;
		cerr << "Error : no tag with name \"" << chemicalElementsDataBasePathTagName << "\" was found." << endl;
		cerr << "Error : no configuration for the chemical elements data base was found" << endl;	
	} else {
		chemicalElementsConfiguration.setChemicalElementsDataBaseRelativePath(getNodeStringValue(chemicalElementsDataBasePathElement,xmlParser::Remove_White_Space));
	}

	DOMNode * chemicalElementsDataBaseSchemaPathElement = getElementByTagName(chemicalElementsDataBaseElement , chemicalElementsDataBaseSchemaPathTagName);

	if( chemicalElementsDataBaseSchemaPathElement== NULL) {
		cerr << "Error : in const pair<string,string> xmlAspicConfigurationParser::_getChemicalElementsDataBaseRelativePath(void) const" << endl;
		cerr << "Error : no tag with name \"" << chemicalElementsDataBaseSchemaPathTagName << "\" was found." << endl;
		cerr << "Error : no configuration for the chemical elements data base was found" << endl;
	} else {
		chemicalElementsConfiguration.setChemicalElementsDataBaseSchemaRelativePath(getNodeStringValue(chemicalElementsDataBaseSchemaPathElement,xmlParser::Remove_White_Space));	
	}

	return chemicalElementsConfiguration;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to find the base URI of the server for ASPIC databases.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
const shellTypesDataBaseConfiguration xmlAspicConfigurationParser::_getShellTypesDataBaseConfiguration(void) const
{	
	const string shellTypesDataBaseTagName = "ShellTypesDataBase";
	const string shellTypesDataBasePathTagName = "ShellTypesDataBasePath";
	const string shellTypesDataBaseSchemaPathTagName = "ShellTypesDataBaseSchemaPath";

	shellTypesDataBaseConfiguration shellTypesConfiguration("","");

	DOMElement * shellTypesDataBaseElement = (DOMElement *)getElementByTagName(shellTypesDataBaseTagName);
	
	if( shellTypesDataBaseElement== NULL) {
		cerr << "Error : in const pair<string,string> xmlAspicConfigurationParser::_getshellTypesDataBaseRelativePath(void) const" << endl;
		cerr << "Error : no tag with name \"" << shellTypesDataBaseTagName << "\" was found." << endl;
		cerr << "Error : no configuration for the chemical elements data base was found" << endl;
		return shellTypesConfiguration;
	}

	DOMNode * shellTypesDataBasePathElement = getElementByTagName(shellTypesDataBaseElement,shellTypesDataBasePathTagName);

	if( shellTypesDataBasePathElement== NULL) {
		cerr << "Error : in const pair<string,string> xmlAspicConfigurationParser::_getshellTypesDataBaseRelativePath(void) const" << endl;
		cerr << "Error : no tag with name \"" << shellTypesDataBasePathTagName << "\" was found." << endl;
		cerr << "Error : no configuration for the chemical elements data base was found" << endl;	
	} else {
		shellTypesConfiguration.setShellTypesDataBaseRelativePath(getNodeStringValue(shellTypesDataBasePathElement,xmlParser::Remove_White_Space));
	}

	DOMNode * shellTypesDataBaseSchemaPathElement = getElementByTagName(shellTypesDataBaseElement , shellTypesDataBaseSchemaPathTagName);

	if( shellTypesDataBaseSchemaPathElement== NULL) {
		cerr << "Error : in const pair<string,string> xmlAspicConfigurationParser::_getshellTypesDataBaseRelativePath(void) const" << endl;
		cerr << "Error : no tag with name \"" << shellTypesDataBaseSchemaPathTagName << "\" was found." << endl;
		cerr << "Error : no configuration for the chemical elements data base was found" << endl;
	} else {
		shellTypesConfiguration.setShellTypesDataBaseSchemaRelativePath(getNodeStringValue(shellTypesDataBaseSchemaPathElement,xmlParser::Remove_White_Space));	
	}

	return shellTypesConfiguration;
}