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
#ifndef _MOLECULE_
#define _MOLECULE_


#include "atom.h"
#include "basisElement.h"
#include "gaussianBasisFunction.h"
#include "gaussianBasisShell.h"
#include <containor.h>


class molecule
{

private :

	/**
	 * An array containing the atoms that compose the molecule.
	 */
	 containor<atom> Atoms;
	
	/**
	 * An array containing the basis data for each element of the molecule.
	 */
	containor<basisElement> BasisElements;

	/**
	 * The charge of the molecule.
	 */
	int Charge;
	
	/**
	 * A description of the molecule.
	 */
	string Description;

	/**
	 * The name of the molecule.
	 */
	string Name;

 
protected:

	/**
	 * Method COPY.
	 */
	void copy(const molecule & m);

public :

	/**
	 * Constructeur par dŽfaut.
	 */
	molecule(void);
	
	/**
	 * Copy constructor.
	 */
	molecule(const molecule & m);
	
	/**
	 * Destructeur.
	 */
	 virtual ~molecule(void);
	
	/**
	 * Method CLEAR.
	 */
	void clear(void);
	
	/**
	 * Method empty.
	 */
	bool empty(void) const;

	/**
   * Method GET for an atom.
	 */
	const atom & getAtom(int item) const;
	
	/**
   * Method GET for an atom.
	 */
	const containor<atom> & getAtoms(void) const;
	
	/**
	 * Method GET for a basis element.
	 */
	const basisElement & getBasis4Atom(int item) const;
	
	/**
	 * Method GET for a basis element.
	 */
	const containor<basisElement> & getBasis4Atoms(void) const;
	
	/**
	 * Method to find the basis name of a given atom.
	 */
	const string & getBasisName4Atom(int atom) const;

	/**
	 * Method GET for the description.
	 */
	const string & getDescription(void) const;

	/**
	 * Méthode pour retrouver l'unité de distance d'un atome.
	 */
	const atom::distanceUnit & getDistanceUnit4Atom(const int & atom) const;

	/**
	 * Method GET for the charge.
	 */
	int getCharge(void) const;
	
	/**
	 * Method GET for a basis function.
	 */
	gaussianBasisFunction getGaussianBasisFunction(const int & atom , const int & shell , const int & function) const;

	/**
	 * Method GET for a basis shell.
	 */
	gaussianBasisShell getGaussianBasisShell(const int & atom , const int & shell) const;

	/**
	 * Method GET for the name.
	 */
	const string & getName(void) const;

	/**
	 * Method GET for the number of atoms.
	 */
	 int getNbrOfAtoms(void) const;
	
	 /**
	  * Method GET for the number of basis function of a given shell.
		*/
	 int getNbrOfBasisFunction4Shell(const int & atom , const int & shell) const;

	 /**
	  * Method GET for the number of basis function of a given shell.
		*/
	 int getNbrOfContractions4Shell(const int & atom , const int & shell) const;

	 /**
	  * Method GET for the number of protons of a given atom.
		*/
	 int getNbrOfProtons4Atom(int atom) const;

	 /**
	 * Method GET for the number of atoms.
	 */
	 int getNbrOfShells4Atom(int item) const;

	/**
	 * Method GET for the number of basis functions.
	 */
	 int getNbrOfBasisFunctions(void) const;

	/**
	 * Method GET for the number of electrons.
	 */
	int getNbrOfElectrons(void) const;
	
	/**
	 * Method GET for the number of protons in the molecule.
	 */
	int getNbrOfProtons(void) const;

	/**
	 * Method GET for the position of a given atom.
	 */
	const dpoint<3> & getPosition4Atom(const int & atom) const;
	
	/**
	 * Method GET for the position of a given expressed in the given unit.
	 */
	const dpoint<3> getPosition4Atom(const int & atom , const atom::distanceUnit & unit) const;
	
	/**
   * Method GET for the unit of in which the position of the center of the atom
	 * is expressed.
   *
	 * @return the distance unit for the position of the center the atom.
   */
	const atom::distanceUnit & getPositionUnit4Atom(const int & atom) const;
	
	/**
	 * Method SET for an atom.
	 */
	void setAtom(const int & item , const atom & a);

	/**
	 * Method SET for an atom.
	 */
	void setAtoms(const containor<atom> & atoms);
	
	/**
	 * Method SET for a basis element.
	 */
	void setBasis4Atom(int item , const basisElement & element);
	
	/**
	 * Method SET for a basis element.
	 */
	void setBasis4Atoms(const containor<basisElement> & basisElements);

	/**
	 * Method SET for the charge of the molecule.
	 */
	void setCharge(int charge);
	
	/**
	 * Method SET for the description.
	 */
	void setDescription(const string & description);
	
	/**
	 * Méthode pour fixer l'unité de distance pour tout les atomes de la molécule.
	 *  Attention cette méthode ne donne pas une valeur à l'unité mais effectue les
	 * conversions nécessaires.
	 *
	 * @param unit l'unité dans laquelle on souhaite exprimer les distances pour la présente molécule.
	 */
	void changeDistanceUnit4Atoms(const atom::distanceUnit & unit) const;

	/**
	 * Method SET for the description.
	 */
	void setName(const string & name);

	/**
	 * Method set for the number of atoms of the molecule.
	 */
	void setNbrOfAtoms(int nbrOfAtoms);

	/**
	 * Method WRITE for the molecule.
	 */
	void write(ostream & out , const string & linePrefix = "") const;


	/**
	 * Method WRITE for the molecule.
	 */
	void writeXML(ostream & out , const string & linePrefix = "") const;

};


/**
 * Extern operator << .
 */
 ostream & operator<<(ostream & out , const molecule & pol);

#endif

