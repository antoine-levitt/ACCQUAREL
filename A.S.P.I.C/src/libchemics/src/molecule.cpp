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
#include "molecule.h"
#include "xmlShellTypesDataBaseInterface.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Default Constructor.
//
// This one constructs an empty molecule.
// The Atoms array is empty (contains 0 elements).
// The BasisElements array is empty (contains 0 elements).
// The Charge is null.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
molecule::molecule(void)
	: Atoms() , BasisElements() , Charge(0) , Description("") , Name("")
{
	;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copy constructor.
//
// This one constructs a copy that is a copy
// of the object of class molcule pass in argument.
//
// @see copy.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
molecule::molecule(const molecule & m)
	: Atoms() , BasisElements() , Charge(0), Description("") , Name("")
{
	copy(m);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Destructor.
//
// This one just clears the current objet.
// @see clear.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
molecule::~molecule(void) 
{
	clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method CLEAR.
//
// This method resets the value of the object.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void molecule::changeDistanceUnit4Atoms(const atom::distanceUnit & unit) const
{
	int i , nbrOfAtoms;
	nbrOfAtoms = getNbrOfAtoms();	
	for(i=0 ; i < nbrOfAtoms ; i++) {
		getAtom(i).changeDistanceUnit(unit);
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method CLEAR.
//
// This method resets the value of the object.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void molecule::clear(void) 
{
	Atoms.clear();
	BasisElements.clear();
	Charge = 0;
	Description = "";
	Name = "";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method COPY.
//
// the purpose of this method is to copy an object of class
// molecule in the current object.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void molecule::copy(const molecule & m)
{
	setNbrOfAtoms(m.getNbrOfAtoms());
	setAtoms(m.getAtoms());
	setBasis4Atoms(m.getBasis4Atoms());
	setCharge(m.getCharge());
	setDescription(m.getDescription());
	setName(m.getName());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to know if the object molecule is empty.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool molecule::empty(void) const
{
	if(getNbrOfAtoms() == 0) {
		return true;
	} else {
		return false;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the atom.
//
// This method returns a const reference to the object
// of class atom stored in the Atom array in position 
// item.
//
// @warning item must be a non negative integer.
//
// @warning item must be strictly inferior than
// the value returned by the method getNbrOfAtoms.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const atom & molecule::getAtom(int item) const
{
	assert(item >= 0);
	assert(item < getNbrOfAtoms());

	return Atoms[item];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the atom's array.
//
// This method returns a const reference to the containor
// of atoms stored in the molecule object.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const containor<atom> & molecule::getAtoms(void) const
{
	return Atoms;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the basisElement.
//
// This method returns a const reference to the object
// of class gaussianBasisElment stored in the basisElement array in position 
// item.
//
// @warning item must be a non negative integer.
//
// @warning item must be strictly inferior than
// the value returned by the method getNbrOfAtoms.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const basisElement & molecule::getBasis4Atom(int item) const
{
	assert(item >= 0);
	assert(item < getNbrOfAtoms());

	return BasisElements[item];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the all the basisElements.
//
// This method returns a const reference to the containor
// of gaussianBasisElment stored in the basisElement array.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const containor<basisElement> & molecule::getBasis4Atoms(void) const
{
	return BasisElements;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to find the basis name for a particular atom.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string & molecule::getBasisName4Atom(int atom) const
{
	assert(atom >= 0);
	assert(atom < getNbrOfAtoms());

	return getBasis4Atom(atom).getBasisName();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the charge.
//
// This method returns the charge of the molecule.
// this value is the difference between the total number
// of protons ans the total number of electrons that are
// present in the molecule.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int molecule::getCharge(void) const
{
	return Charge;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the description of the molecule.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string & molecule::getDescription(void) const
{
	return Description;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour retrouver l'unité de distance dans laquelle est exprimé la position d'un atome.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const atom::distanceUnit & molecule::getDistanceUnit4Atom(const int & atom) const
{
	assert(atom >= 0);
	assert(atom < getNbrOfAtoms() );

	return getAtom(atom).getDistanceUnit();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode qui construit la fonction de base pour le triplet atom - couche - fonction de base.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
gaussianBasisFunction molecule::getGaussianBasisFunction(const int & atom , const int & shell , const int & function) const
{
	assert(atom >= 0);
	assert(atom < getNbrOfAtoms() );

	assert(shell >=0 );
	assert(shell < getNbrOfShells4Atom(atom));

	assert(function >= 0);
	assert(function < getNbrOfBasisFunction4Shell(atom , shell));

	gaussianBasisFunction basisFunction;
	
	
		//* DEPRECATED ....
		//Now xe directly ask the right unit to the atom molecule class.
		
		dpoint<3> position;
		// Bon les fonctions de base gaussienne on des coefficients 
	// en (unité atomiques)^(-2). Il faut donc convertir les positions
	// sinon les calculs effectuées sur les fonctions de base sont faux.
	if(getPositionUnit4Atom(atom) != atom::ATOMIC_UNIT) {
			
		// ici je mets un assert qui ne sert pas sauf que si on veut rajouter une autre unité il faudra faire autrement.
		// et cela aura au moins le mérite d'arreter le programme à la génèse de l'erreur.
		assert(getPositionUnit4Atom(atom) == atom::ANGSTROEM);
		
		position = atom::angstroem2AtomicUnitConverter(getPosition4Atom(atom));	
	} else {
		position = getPosition4Atom(atom);
	}

	// une fois les conversion effectuée il ne reste plus qu'à recopier la valeur 
	// de la position.
	basisFunction.setCenter(position);
	
	//*/
	//basisFunction.setCenter(getPosition4Atom(atom,atom::ATOMIC_UNIT));
	basisFunction.setContractions(getBasis4Atom(atom).getContractions4Shell(shell));
	basisFunction.setMonomeDegree(getBasis4Atom(atom).getMonomeDegree(shell,function));

	return basisFunction;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode qui construit la fonction de base pour la paire atome - couche.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
gaussianBasisShell molecule::getGaussianBasisShell(const int & atom , const int & shell) const
{
	assert(atom >= 0);
	assert(atom < getNbrOfAtoms() );

	assert(shell >=0 );
	assert(shell < getNbrOfShells4Atom(atom));

	gaussianBasisShell basisShell;
	
	
	
	
	/* DEPRECATED ....
		Now xe directly ask the right unit to the atom molecule class.
	
	 dpoint<3> position;
	// Bon les fonctions de base gaussienne on des coefficients 
	// en (unité atomiques)^(-2). Il faut donc convertir les positions
	// sinon les calculs effectuées sur les fonctions de base sont faux.
	if(getPositionUnit4Atom(atom) != atom::ATOMIC_UNIT) {
		
		// ici je mets un assert qui ne sert pas sauf que si on veut rajouter une autre unité il faudra faire autrement.
		// et cela aura au moins le mérite d'arreter le programme à la génèse de l'erreur.
		assert(getPositionUnit4Atom(atom) == atom::ANGSTROEM);
		
		position = atom::angstroem2AtomicUnitConverter(getPosition4Atom(atom));	
	} else {
		position = getPosition4Atom(atom);
	}

	// une fois les conversion effectuée il ne reste plus qu'à recopier la valeur 
	// de la position.
	basisShell.setCenter(position);
	*/
	
	basisShell.setCenter(getPosition4Atom(atom , atom::ATOMIC_UNIT));
	
	// On met ensuite les contractions
	basisShell.setContractions(getBasis4Atom(atom).getContractions4Shell(shell));
	
	// On met ensuite le type de couche
	basisShell.setShellType(getBasis4Atom(atom).getShellType(shell));

	return basisShell;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the name of the molecule.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string & molecule::getName(void) const
{
	return Name;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the number of atoms in the molecule.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int molecule::getNbrOfAtoms(void) const
{
	return Atoms.getSizes();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the number of basis functions.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int molecule::getNbrOfBasisFunctions(void) const
{
	int i , j , nbrOfAtoms , nbrOfShells , nbrOfBasisFunctions =0;

	nbrOfAtoms = getNbrOfAtoms();

	for(i=0 ; i < nbrOfAtoms ; i++) {
		nbrOfShells = getBasis4Atom(i).getNbrOfShells();
		for(j=0 ; j < nbrOfShells ; j++) {
			nbrOfBasisFunctions += getNbrOfBasisFunction4Shell(i,j);
		}
	}

	return nbrOfBasisFunctions;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the number of basis functions of a given shell.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int molecule::getNbrOfBasisFunction4Shell(const int & atom , const  int & shell) const
{
	assert(atom >= 0);
	assert(atom < getNbrOfAtoms());

	assert(shell >= 0);
	assert(shell < getNbrOfShells4Atom(atom));

	return getBasis4Atom(atom).getNbrOfBasisFunctions4Shell(shell);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the number of basis functions of a given shell.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int molecule::getNbrOfContractions4Shell(const int & atom , const int & shell) const
{
	assert(atom >= 0);
	assert(atom < getNbrOfAtoms());

	assert(shell >= 0);
	assert(shell < getNbrOfShells4Atom(atom));

	return getBasis4Atom(atom).getNbrOfContractions4Shell(shell);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the number of electrons in the molecule.
//
// the number of electrons is computed using the facts that
// this number is the total number of protons soustracted with
// the charge of a the molecule.
//
// @warning if the charge is equal to the total number of protons
// the molecule does not have any electrons. This is possible but
// generates a warning.
//
// @warning if the charge is strictly higher than the total number
// of protons, the molecule would have a negative number of electrons
// This is impossible and generates an error.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int molecule::getNbrOfElectrons(void) const
{
	int nbrOfElectrons;

	nbrOfElectrons = getNbrOfProtons() - getCharge();

	if(nbrOfElectrons == 0) {
		cerr << "Warning : in int molecule::getNbrOfElectrons(void) const" << endl;
		cerr << "Warning : no electrons are in the molecule." << endl;
		cerr << "Warning : make sure this is what you want !" << endl;
	}	

	if(nbrOfElectrons < 0) {
		cerr << "Error : in int molecule::getNbrOfElectrons(void) const" << endl;
		cerr << "Error : a negative number of electrons was found." << endl;
		cerr << "Error : this is a fatal error." << endl;
		exit(1);
	}

	return nbrOfElectrons;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the total number of protons in the molecule.
//
// the number of protons is the sum on all atoms of the number 
// of protons of the atom.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int molecule::getNbrOfProtons(void) const
{
	int atom , nbrOfAtoms , nbrOfProtons;

	nbrOfProtons = 0;
	nbrOfAtoms = getNbrOfAtoms();

	for(atom = 0 ; atom < nbrOfAtoms ; atom++) {
		nbrOfProtons += getAtom(atom).getNbrOfProtons();
	}

	return nbrOfProtons;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the number of shell for a given atom.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int molecule::getNbrOfProtons4Atom(int atom) const
{
	assert(atom >= 0);
	assert(atom < getNbrOfAtoms());

	return getAtom(atom).getNbrOfProtons();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the number of shell for a given atom.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int molecule::getNbrOfShells4Atom(int atom) const
{
	assert(atom >= 0);
	assert(atom < getNbrOfAtoms());

	return getBasis4Atom(atom).getNbrOfShells();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the position of a given atom.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const dpoint<3> &  molecule::getPosition4Atom(const int & atom) const
{
	assert(atom >= 0);
	assert(atom < getNbrOfAtoms());

	return getAtom(atom).getPosition();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the position of a given atom.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const dpoint<3>  molecule::getPosition4Atom(const int & atom , const atom::distanceUnit & unit) const
{
	assert(atom >= 0);
	assert(atom < getNbrOfAtoms());
	
	return getAtom(atom).getPosition(unit);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the position of a given atom.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const atom::distanceUnit &  molecule::getPositionUnit4Atom(const int & atom) const
{
	assert(atom >= 0);
	assert(atom < getNbrOfAtoms());

	return getAtom(atom).getDistanceUnit();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method SET for an atom of the molecule.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void molecule::setAtom(const int & item , const atom &  a)
{
	assert( item >= 0);
	assert( item < getNbrOfAtoms());

	Atoms[item] = a;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method SET for all atoms of the molecule.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void molecule::setAtoms(const containor<atom> &  atoms)
{
	assert(getNbrOfAtoms() == atoms.getSizes());
	Atoms = atoms;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method SET for a basis of an atom of the molecule.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void molecule::setBasis4Atom(int item , const basisElement & basisElement)
{
	assert( item >= 0);
	assert( item < getNbrOfAtoms());
	
	BasisElements[item] = basisElement;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method SET for all basis of an atom of the molecule.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void molecule::setBasis4Atoms(const containor<basisElement> & basisElements)
{
	assert(basisElements.getSizes() == getNbrOfAtoms());
	
	BasisElements = basisElements;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method SET for all basis of an atom of the molecule.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void molecule::setCharge(int charge)
{
	Charge = charge;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method SET for all basis of an atom of the molecule.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void molecule::setDescription(const string & description)
{
	Description = description;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method SET for all basis of an atom of the molecule.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void molecule::setName(const string & name)
{
	Name = name;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method SET for the number of atoms.
//
// This method allocates space for the atoms and the basis
// elements that will be stored in the object.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void molecule::setNbrOfAtoms(int nbrOfAtoms) 
{
	assert(nbrOfAtoms > 0);

	Atoms.setSizes(nbrOfAtoms);
	BasisElements.setSizes(nbrOfAtoms);
}

void molecule::write(ostream & outStream , const string & linePrefix) const
{
	string atomLinePrefix = "\t" + linePrefix;

	outStream << linePrefix << " - Name        : " << getName() << "." << endl;
	outStream << linePrefix << " - Description : " << getDescription() << endl; 
	outStream << "---------------------------------------------------------------" << endl; 	
	
	for(int i=0 ; i < getNbrOfAtoms() ; i++) {
		outStream << linePrefix << " - Atom " << i << "." << endl;
		getAtom(i).write(outStream,atomLinePrefix);
		outStream << atomLinePrefix << " - Basis             : \"" << getBasisName4Atom(i) << "\"." << endl; 
	outStream << "---------------------------------------------------------------" << endl; 	
	}
	outStream << linePrefix << " - Number of Protons        : " << getNbrOfProtons() << "." << endl;
	outStream << linePrefix << " - Number of electrons      : " << getNbrOfElectrons() << "." << endl;
	outStream << linePrefix << " - Charge                   : " << getCharge() << "." << endl;
	outStream << linePrefix << " - Number of Basis Function : " << getNbrOfBasisFunctions() << "." << endl;
}


void molecule::writeXML(ostream & outStream , const string & linePrefix) const
{
	string atomLinePrefix = linePrefix + "\t\t\t";

	outStream << "<Molecule>" << endl;
	outStream << "\t<Name>"<< getName() << "</Name>" << endl;
	outStream << "\t<Description>" << getDescription() << "</Description>" << endl; 
	outStream << linePrefix << "\t<Charge>" << getCharge() << "</Charge>" << endl;
	outStream << linePrefix << "\t<ElementsList>" << endl;
	for(int i=0 ; i < getNbrOfAtoms() ; i++) {
		outStream << linePrefix << "\t\t<Element>" << endl;
		getAtom(i).writeXML(outStream,atomLinePrefix);
		outStream << atomLinePrefix << "<BasisKey>" << getBasisName4Atom(i) << "</BasisKey>" << endl; 
		outStream << linePrefix << "\t\t</Element>" << endl;
	}
	outStream << linePrefix << "\t</ElementsList>" << endl;
	outStream << linePrefix << "</Molecule>" << endl;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Extern Operator << for the molecule objects.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ostream & operator<<(ostream & outStream , const molecule & mol)
{
	mol.write(outStream);
	return outStream;
}

