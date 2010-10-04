/* 
 * Sample program of the A.S.P.I.C plate forme :
 * This program generates HTML files from the data bases
 * in ASPIC such that users can refer to an web page instead 
 * of an XML file to view the informations.
 *
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
#include <aspicConfiguration.h>
#include <chemicalElementsDataBaseInterface.h>
#include <xmlBasisDataBaseInterface.h>
#include <iostream>
#include <fstream>
#include <xmlShellTypesDataBaseInterface.h>
#include <string4Polynomes.h>
#include <xmlBasisKeyListParser.h>
using namespace std;

/**
 * Fonction qui va créer le fichier HTML pour la base de données des types de couches.
 */
void basisPrinter(void) 
{
	bool newRow;
	int nbrOfBasis , i , nbrOfElements , j  , nbrOfContractions , k , nbrOfShells , l;
	basisElement basisElement;
	containor<string> basisList;
	
	
	string outFileName = aspicConfiguration::getAspicRoot();
	outFileName += "/documentation/gaussianBasis.html";
	
	ofstream out;
	out.open(outFileName.c_str());

	if(out.is_open() == false) {
		cerr << "Error : in void basisPrint(void)" << endl;
		cerr << "Error : unable to open file " << outFileName << " for writing." << endl;
		cerr << "Error : aborting ..." << endl;
		return;
	}



	basisList = xmlBasisKeyListParser::getBasisKeyList();
	nbrOfBasis = basisList.getSizes();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Les en têtes du document html.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "<html>" << endl;
	out << "<head>" << endl;
	out << "<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">" << endl;
	out << "<title>The Basis Data Base of A.S.P.I.C</title>" << endl;
	out << "<style type=\"text/css\">"<< endl;
	out << "h1 {color: red}" << endl;
	out << "h2 {color: blue}" << endl;
	out << "h3 {color: green}" << endl;
	out << "th.bg {background-color: #F0F8FF; text-align: left}" << endl;
	out << "</style>" << endl;

	out << "</head>" << endl;

	out << "<body>" << endl;
	out << "<table style=\"width: 100%\">" << endl;
	out << "<tr>" << endl;
	out << "<td>The <a href=\"shellTypes.html\">shell types</a> Data Base</td>" << endl;
	out << "<td style=\"text-align: center\"><a href=\"index.html\">Home</a></td>" << endl;
	out << "<td style=\"text-align: right\">The <a href=\"chemicalElements.html\">Chemical Elements</a> Data Base</td>" << endl;
	out << "</tr>" << endl;
	out << "</table>" << endl;
	out << "<div style=\"text-align: center;\">" << endl;
	out << "<h1><big style=\"font-weight: bold;\">The Basis Data Base of A.S.P.I.C.</big></h1>" << endl;
	out << "<h1>(" << nbrOfBasis << " Basis)</h1></div>" << endl;
	out << "<br>" << endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Tableau du début avec tout les symobles présents dans la base de données.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "<a name=\"BasisList\"><h2>The Basis List</h2></a>" << endl;
	newRow = true;
	out << "<table border=\"1\" cellpadding=\"10\">" << endl;

	for(i=0 ; i<nbrOfBasis ; i++) {
		
		if(newRow){
			out << "<tr>" << endl;
			newRow = false;
		}

		out << "<td><a href=\"#Basis-"<< i << "\">" << basisList[i] << "</a></td>" << endl;
		
		if(i%7==6) {
			out << "</tr>" << endl;
			newRow = true;
		}
		
	}
	out << "</table>" << endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Tableau avec toute les bases
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(i=0 ; i < nbrOfBasis ; i++) {
		
		xmlBasisDataBaseInterface::connect(basisList[i]);
		nbrOfElements = xmlBasisDataBaseInterface::getNbrOfBasisElements();

		out << "<a name=\"Basis-" << i << "\"><h2>Basis \"" << basisList[i] << "\"" << endl;
		out << "("<< (i+1) <<"/" << nbrOfBasis << ")"<< endl;
		out << " contains " << nbrOfElements << " Chemical Element(s)." << endl;
		out << "</h2></a>"<< endl;
	
		out << "This file was generated by the aspicDataBasePrinter program from the <a href=\"../data/chemics/gaussianBasis/" << basisList[i] << ".xml\">gaussian basis</a> Data Base. " << endl;
		out << "<b>Warning :</b> this link is points to the unformated XML file that contains the data base." << endl;
		out << " Make sure your browser is compatible or use the save as option." << endl;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Table avec tout les éléments de la base.
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		out << "<a name=\"Elements-" << i <<"\"><h3>The Element List</h3></a>" << endl;
		newRow = true;
		out << "<table border=\"1\" cellpadding=\"10\">" << endl;
		for(j=0 ; j<nbrOfElements ; j++) {
			basisElement = xmlBasisDataBaseInterface::getBasisElement(j);
			if(newRow){
				out << "<tr>" << endl;
				newRow = false;
			}
			out << "<td><a href=\"#Basis-"<< i << "-Element-" << j <<"\">" << basisElement.getBasisElementKey() << "</a></td>" << endl;
			if(i%7==6) {
				out << "</tr>" << endl;
				newRow = true;
			}

		}
		out << "</table>" << endl;
		out << "Back to the <a href=\"#BasisList\">basis list</a>." << endl;
		
	
		for(j=0 ; j<nbrOfElements ; j++) {
			basisElement = xmlBasisDataBaseInterface::getBasisElement(j);
			nbrOfShells = basisElement.getNbrOfShells();
			out << "<a name=\"Basis-"<< i << "-Element-" << j <<"\"><h3>" << basisElement.getBasisElementKey() << "</h3></a>" << endl;
			
			out << "<table border=\"1\">" << endl;
			
			for(k=0; k < nbrOfShells; k++) {
				out << "<tr><th class=\"bg\" align=\"center\"> Shell " << (k+1) << "/" << nbrOfShells << "</td>" << endl;
				out << "<td>" << endl;		
				out << "<table border=\"1\">" << endl;		
				out << "<tr>" << endl;
				out << "<th class=\"bg\">Shell Type</th>" << endl;
				out << "<td align=\"center\">" << basisElement.getShellTypeKey(k) << "</td>" << endl;
				out << "</tr>" << endl;
				out <<"<tr>" << endl;
				out << "<th class=\"bg\">Contractions</th>" << endl;
				out << "<td>" << endl;
				out << "<table>" << endl;
				out << "<tr><th>Coefficient</th><th>Exponent</th></tr>" << endl;
				nbrOfContractions = basisElement.getContractions4Shell(k).getNbrOfContractions();
				for(l=0 ; l < nbrOfContractions ; l++) {
				out << "<tr>" << endl;
				out << "<td>" << basisElement.getContractions4Shell(k).getCoefficient(l) <<"</td>" << endl; 
				out << "<td>" << basisElement.getContractions4Shell(k).getExponent(l) << "</td></tr>" << endl;
				}
				out << "</table>" << endl;
				out << "</td>"		<< endl;
				out <<"</tr>" << endl;
				out << "</table>" << endl;		
				out << "</td></tr>" << endl;
			}
			out << "</table>" << endl;
			out << "Back to the <a href=\"#BasisList\">basis list</a>.<br>" << endl;
			out << "Back to the <a href=\"#Elements-"<< i << "\">element list</a>." << endl;
			

		}

	}

}

/**
 * Fonction qui va créer le fichier HTML pour la base de données des types de couches.
 */
void shellTypesPrinter(void) {

	bool newRow;
	int nbrOfShellTypes , i , nbrOfBasisFunctions , j ;
	shellType basisShell;
	
	
	
	string outFileName = aspicConfiguration::getAspicRoot();
	outFileName += "/documentation/shellTypes.html";
	
	ofstream out;
	out.open(outFileName.c_str());

	if(out.is_open() == false) {
		cerr << "Error : in void shellTypesPrint(void)" << endl;
		cerr << "Error : unable to open file " << outFileName << " for writing." << endl;
		cerr << "Error : aborting ..." << endl;
		return;
	}
	
	

	nbrOfShellTypes = xmlShellTypesDataBaseInterface::getNbrOfShellTypes();
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Les en têtes du document html.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "<html>" << endl;
	out << "<head>" << endl;
	out << "<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">" << endl;
	out << "<title>The Shell Types Data Base of A.S.P.I.C</title>" << endl;
	out << "<style type=\"text/css\">"<< endl;
	out << "h1 {color: red}" << endl;
	out << "h2 {color: blue}" << endl;
	out << "th.bg {background-color: #F0F8FF; text-align: left}" << endl;
	out << "</style>" << endl;

	out << "</head>" << endl;

	out << "<body>" << endl;
	out << "<table style=\"width: 100%\">" << endl;
	out << "<tr>" << endl;
	out << "<td>The <a href=\"gaussianBasis.html\">Gaussian Basis</a> Data Base</td>" << endl;
	out << "<td style=\"text-align: center\"><a href=\"index.html\">Home</a></td>" << endl;
	out << "<td style=\"text-align: right\">The <a href=\"chemicalElements.html\">Chemical Elements</a> Data Base</td>" << endl;
	out << "</tr>" << endl;
	out << "</table>" << endl;
	out << "<div style=\"text-align: center;\">" << endl;
	out << "<h1><big style=\"font-weight: bold;\">The Shell Types Data Base of A.S.P.I.C.</big></h1>" << endl;
	out << "<h1>(" << nbrOfShellTypes << " Types)</h1></div>" << endl;
	out << "<br>" << endl;
	out << "This file was generated by the aspicDataBasePrinter program from the <a href=\"../data/chemics/shellTypes.xml\">Shell Types</a> Data Base. " << endl;
	out << "<b>Warning :</b> this link is points to the unformated XML file that contains the data base." << endl;
	out << " Make sure your browser is compatible or use the save as option." << endl;
	

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Tableau du début avec tout les symobles présents dans la base de données.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "<a name=\"List\"><h2>The Shell Type List</h2></a>" << endl;
	newRow = true;
	out << "<table border=\"1\" cellpadding=\"10\">" << endl;

	for(i=0 ; i<nbrOfShellTypes ; i++) {
		basisShell = xmlShellTypesDataBaseInterface::getShellType(i);
		
		
		if(newRow){
			out << "<tr>" << endl;
			newRow = false;
		}

		out << "<td><a href=\"#Element-"<< i << "\">" << basisShell.getShellTypeKey() << "</a></td>" << endl;
		
		if(i%7==6) {
			out << "</tr>" << endl;
			newRow = true;
		}
		
	}
	out << "</table>" << endl;



	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Pour chaque élément on écrit un tat de chose ...
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(i=0 ; i< nbrOfShellTypes ; i++) {
	
		basisShell = xmlShellTypesDataBaseInterface::getShellType(i);
		
		out << "<a name=\"Element-" << i << "\"><h2>Element \"" << basisShell.getShellTypeKey() <<"\" ("<< (i+1) <<"/" << nbrOfShellTypes << ")</h2></a>"<< endl;
		out << "<table border=\"1\">" << endl;
		
		out << "<tr>" << endl;
		out << "<th class=\"bg\"> Key : </th>" << endl;
		out << "<td> \"" << basisShell.getShellTypeKey() << "\"</td>" << endl;
		out << "</tr>" << endl;

		out << "<tr>" << endl;
		out << "<th class=\"bg\"> Name : </th>" << endl;
		out << "<td>\"" << basisShell.getShellTypeName() << "\"</td>" << endl;
		out << "</tr>" << endl;

		out << "<tr>" << endl;
		out << "<th class=\"bg\"> Monome Degrees : </th>" << endl;
		out << "<td>" << endl;
		
		out<< "<table border=\"1\">" << endl;
		out <<"<tr><th align=\"center\">x</th><th align=\"center\">y</th><th align=\"center\">z</th></tr>" << endl;
		nbrOfBasisFunctions = basisShell.getNbrOfBasisFunctions();
		for( j=0 ; j < nbrOfBasisFunctions ; j++) {
			out << "<tr>" << endl;
			out << "<td align=\"center\">" << basisShell.getMonomeDegree(j)[0] << "</td>" << endl;
			out << "<td align=\"center\">" << basisShell.getMonomeDegree(j)[1] << "</td>" << endl;
			out << "<td align=\"center\">" << basisShell.getMonomeDegree(j)[2] << "</td>" << endl;
			out << "</tr>" << endl;

		}
		out<<	"</table>" << endl;
		
		out << "</td>" << endl;
		out << "</tr>" << endl;

		out << "</table>" << endl;
		out << "return to the <a href=\"#List\">list of all the shell types</a>.<br>" << endl;
		
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// La fin du document HTML.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "</body>" << endl;
	out << "</html>" << endl;

}



/**
 * Fonction qui va créer le fichier HTML pour la base de données des éléments chimiques.
 */
void chemicalElementsPrinter(void) 
{

	bool newRow;
	int nbrOfChemicalElements , i , nbrOfIsotopes , j ;
	chemicalElement element;
	
	
	
	string outFileName = aspicConfiguration::getAspicRoot();
	outFileName += "/documentation/chemicalElements.html";
	
	ofstream out;
	out.open(outFileName.c_str());

	if(out.is_open() == false) {
		cerr << "Error : in void chemicalElementsPrint(void)" << endl;
		cerr << "Error : unable to open file " << outFileName << " for writing." << endl;
		cerr << "Error : aborting ..." << endl;
		return;
	}
	
	nbrOfChemicalElements = xmlChemicalElementsDataBaseInterface::getNbrOfChemicalElements();
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Les en têtes du document html.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "<html>" << endl;
	out << "<head>" << endl;
	out << "<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">" << endl;
	out << "<title>The Chemical Elements Data Base of A.S.P.I.C</title>" << endl;
	out << "<style type=\"text/css\">"<< endl;
	out << "h1 {color: red}" << endl;
	out << "h2 {color: blue}" << endl;
	out << "th.bg {background-color: #F0F8FF; text-align: left}" << endl;
	out << "</style>" << endl;

	out << "</head>" << endl;

	out << "<body>" << endl;
	out << "<table style=\"width: 100%\">" << endl;
	out << "<tr>" << endl;
	out << "<td>The <a href=\"gaussianBasis.html\">Gaussian Basis</a> Data Base</td>" << endl;
	out << "<td style=\"text-align: center\"><a href=\"index.html\">Home</a></td>" << endl;
	out << "<td style=\"text-align: right\">The <a href=\"shellTypes.html\">Shell Types</a> Data Base</td>" << endl;
	out << "</tr>" << endl;
	out << "</table>" << endl;
	out << "<div style=\"text-align: center;\">" << endl;
	out << "<h1><big style=\"font-weight: bold;\">The Chemical Elements Data Base of A.S.P.I.C.</big></h1>" << endl;
	out << "<h1>(" << nbrOfChemicalElements << " Elements)</h1></div>" << endl;
	out << "<br>" << endl;
	out << "This file was generated by the aspicDataBasePrinter program from the <a href=\"../data/chemics/chemicalElements.xml\">Chemical Elements Data Base</a>. " << endl;
	out << "<b>Warning :</b> this link is points to the unformated XML file that contains the data base." << endl;
	out << " Make sure your browser is compatible or use the save as option." << endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Tableau du début avec tout les symobles présents dans la base de données.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "<a name=\"List\"><h2>The Elements List</h2></a>" << endl;
	newRow = true;
	out << "<table border=\"1\" cellpadding=\"10\">" << endl;

	for(i=0 ; i<nbrOfChemicalElements ; i++) {
		element = xmlChemicalElementsDataBaseInterface::getChemicalElement(i);
		
		
		if(newRow){
			out << "<tr>" << endl;
			newRow = false;
		}

		out << "<td><a href=\"#Element-"<< i << "\">" << element.getChemicalElementKey() << "</a></td>" << endl;
		
		if(i%7==6) {
			out << "</tr>" << endl;
			newRow = true;
		}
		
	}
	out << "</table>" << endl;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Pour chaque élément on écrit un tat de chose ...
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(i=0 ; i<nbrOfChemicalElements ; i++) {
		element = xmlChemicalElementsDataBaseInterface::getChemicalElement(i);
		out << "<a name=\"Element-" << i << "\"><h2>Element \"" << element.getSymbol() <<"\" ("<< (i+1) <<"/" << nbrOfChemicalElements << ")</h2></a>"<< endl;
		out << "<table border=\"1\">" << endl;
		
		out << "<tr>" << endl;
		out << "<th class=\"bg\"> Key : </th>" << endl;
		out << "<td> \"" << element.getChemicalElementKey() << "\"</td>" << endl;
		out << "</tr>" << endl;

		out << "<tr>" << endl;
		out << "<th class=\"bg\"> Symbol : </th>" << endl;
		out << "<td>\"" << element.getSymbol() << "\"</td>" << endl;
		out << "</tr>" << endl;

		out << "<tr>" << endl;
		out << "<th class=\"bg\"> Name : </th>" << endl;
		out << "<td>\"" << element.getName() << "\"</td>" << endl;
		out << "</tr>" << endl;

		out << "<tr>" << endl;
		out << "<th class=\"bg\"> Nbr of protons : </th>" << endl;
		out << "<td>" << element.getNbrOfProtons() << "</td>" << endl;
		out << "</tr>" << endl;

		out << "<tr>" << endl;
		out << "<th class=\"bg\"> Isotopes : </th>" << endl;
		out << "<td>" << endl;
		
		out<< "<table border=\"1\">" << endl;
		out <<"<tr><th align=\"left\">Nbr of neutrons</th><th align=\"left\">Probability</th></tr>" << endl;
		nbrOfIsotopes = element.getNbrOfIsotopes();
		for( j=0 ; j < nbrOfIsotopes ; j++) {
			out << "<tr>" << endl;
			out << "<td align=\"center\">" << element.getNbrOfNeutrons4Isotope(j) << "</td>" << endl;
			out << "<td align=\"center\">" << element.getProbability4Isotope(j) << "</td>" << endl;
			out << "</tr>" << endl;

		}
		out<<	"</table>" << endl;
		
		out << "</td>" << endl;
		out << "</tr>" << endl;

		out << "</table>" << endl;
		out << "return to the <a href=\"#List\">list of all the chemical elements</a>.<br>" << endl;
		
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// La fin du document HTML.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "</body>" << endl;
	out << "</html>" << endl;
}


int main(int argc , char ** argv) 
{
	basisPrinter();
	chemicalElementsPrinter();
	shellTypesPrinter();
	return 0;
}

