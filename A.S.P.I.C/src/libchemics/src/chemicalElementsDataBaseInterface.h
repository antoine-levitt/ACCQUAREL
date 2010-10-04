/* 
* The chemics library of A.S.P.I.C. 
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
#ifndef _CHEMICAL_ELEMENTS_DATABASE_INTERFACE_
#define _CHEMICAL_ELEMENTS_DATABASE_INTERFACE_

#include <xmlDataBaseInterface.h>
#include "chemicalElement.h"

/**
 * The purpose of this class retrieve chemical elements data
 * in a xml file.
 *
 * The files that are manage by this class can be found in 
 * $(ASPIC_ROOT)/chemicsDataBase/chemicalElements. This file
 * contains data about chemical element, and the searhing key is 
 * usually the symbol of the element.
 */
class xmlChemicalElementsDataBaseInterface : public xmlDataBaseInterface
{
private:
	/**
	 * A pointer to the only object chemical data base interface that will be instanciated. 
	 */
	static xmlChemicalElementsDataBaseInterface * ChemicalElementsDataBase;
	
	/**
	 * The name of the tag that separetes chemical elements in the data base.
	 */
	static const string ChemicalElementKeyAttributeName;
	
	/**
	 * The tag name containing the name of the chemical element.
	 * In the example below this variable should contains "Name".
	 */ 
	static const string ChemicalElementNameTagName;
	
	/**
	 * The tag name containing the number of protons of the 
	 * chemical element.
	 * In the example below this variable should contains "NbrOfProtons".
	 */ 
	static const string ChemicalElementNbrOfProtonsTagName;
	
	/**
	 * The tag name containing the symbol of the chemical element.
	 * In the example below this variable should contains "Symbol".
	 */ 
	static const string ChemicalElementSymbolTagName;
	
	/**
	 * The name of the tag that separetes chemical elements in the data base.
	 */
	static const string ChemicalElementTagName;

	/**
	 * Le nom du tag qui contient les isotopes.
	 */
	static const string ChemicalElementIsotopeTagName;

	/**
	 * Le nom du tag qui contient le nombre de neutrons d'un isotope.
	 */
	static const string ChemicalElementIsotopeNbrOfNeutronsTagName;

	/**
	 * Le nom du tag qui contient la probabilité d'un isotope.
	 */
	static const string ChemicalElementIsotopeProbabilityTagName;

protected:
	
	/**
	 * Method GET that constructs the URI of the document to parse.
	 *
	 * The purpose is to find the file ChemicalElementsBaseName.xml. When working
	 * localy the URI of the file is : $(ASPIC_ROOT)/chemicsDataBase/chemicalElements/chemicalElements.xml.
	 * For the moment there is only this one file, but this method will be convinient 
	 * if one wants to make A.S.P.I.C. more customisable.
	 *
	 * @return the location of the file that contains the data to be parsed.
	 */
	virtual string getDocumentURI(void) const;
	
	/**
	 * Method GET that constructs the URI of the schema use to validate
	 * the document to be parsed.
	 *
	 * The purpose is to find the file ChemicalElements.xsd. The file 
	 * is genraly located in $(ASPIC_ROOT)/chemicsDataBase/chemicalElements/chemicalElements.xsd.
	 *
	 * @return the location of the validating schema.
	 */
	virtual string getSchemaURI(void) const;

	/**
	 * Method GET for a chemical Element.
	 */
	chemicalElement getElement(DOMElement * rootElement) const;
	
	/**
		* The constructor.
	 */
	xmlChemicalElementsDataBaseInterface(void);
	
	/**
	 * The destructor.
	 */
	virtual ~xmlChemicalElementsDataBaseInterface(void);
	
public:

	/**
	 * Method to connect to the chemical elements data base.
	 */		
	static void connect(void);

	/**
	 * Method to connect to the chemical elements data base.
	 */		
	static bool connected(void);
	
	/**
	 * Method to disconnect from the chemical elements data base. 
	 */
	static void disconnect(void);
	
	/**
	 * Method GET.
	 *
	 * @param elementKey the key to search in the data base. Usually this
	 * is the symbol of the element we want to know about.
	 *
	 * @param element a chemicalElement object that will recieve the information
	 * written in the data base.
	 */
	static chemicalElement getChemicalElement(const string & elementKey);

	/**
	 * Method GET.
	 */
	static chemicalElement getChemicalElement(int item);

	/**
	 * Method GET for the number of chemical elements in the data base.
	 */
	static int getNbrOfChemicalElements(void);
};



#endif

