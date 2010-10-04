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
#ifndef _MOLECULE_MAPPER_
#define _MOLECULE_MAPPER_

#include <containor.h>
#include "contractions.h"
#include "molecule.h"
#include "gaussianBasisFunction.h"
#include "gaussianBasisShell.h"
#include "shellType.h"

class basisFunctionPosition
{
public:
	/**
	 * The number of the atom the basis function is 
	 * attached to.
	 */
	int Atom;

	/**
	 * The number of the atom the basis function is 
	 * attached to.
	 */
	int Shell;

	/**
	 * The number of the basis function.
	 */
	int BasisFunction;
};


class moleculeMapper {

private:
		
	/**
	 * Molecule that this numerotation is attach to.
	 */
	const molecule * Molecule;

	/**
	 * This array is used for the basis functions numertation.
	 */
	containor<int> BasisFunctions4Shells;

	/**
	 * This is the number of atoms that are in the molecule.
	 */
	int NbrOfAtoms;

	/**
	 * This array is use to get the shell numerotation.
	 */
	containor<int> Shells4Atoms;

protected:

	/**
	 * Method that computes the numerotation tables.
	 */
	void computeBasisFunctionsTable(void);

	/**
	 * Method that computes the numerotation tables.
	 */
	void computeShellsTable(void);

	/**
	 * Method that computes the numerotation tables.
	 */
	void computeTables(void);

	/**
	 * Method for copy.
	 */
	void copy(const moleculeMapper & map);

	/**
	 * Method SET for the molecule the map is associated to.
	 */
	void setMolecule(const molecule * mol);

public:

	/**
	 * Default constructor.
	 */
	moleculeMapper(void);

	/**
	 * Destructor.
	 */
	~moleculeMapper(void);

	/**
	 * Method to attach a particular molecule to the map.
	 */
	void attach(const molecule & mol);

	/**
	 * Method clear.
	 */
	void clear(void);
	
	/**
	 * Method to know if the map is empty.
	 */
	bool empty(void) const;

	/**
	 * Method GET for the basis functions numerotation table.
	 */
	const containor<int> & getBasisFunctions4Shells(void) const;

	/**
	 * Method GET for the global number of a shell.
	 */
	int getGlobalBasisFunctionNumber(int atom , int shell , int basisFunction) const;
	
	/**
	 * Method GET for the global number of a shell.
	 */
	int getGlobalShellNumber(int atom , int shell) const;
		
	/**
	 * Method GET for the basis function with a global number.
	 */
	gaussianBasisFunction getBasisFunction4Global(int basisFunction) const;

	/**
	 * Method GET for the basis function with a global number.
	 */
	const dpoint<3> & getBasisFunctionCenter(int basisFunction) const;

	/**
	 * Method GET for the basis function with a global number.
	 */
	const dpoint<3>  & getShellCenter(int basisFunction) const;

	/**
	 * Method GET to find the basis function position from a global number of 
	 * basis function.
	 */
	basisFunctionPosition getBasisFunctionPosition4Global(int basisFunctionNumber) const;

	/**
   * Method pour connaitre le nombre de contractions pour une couche.
	 */
	int getNbrOfContractions4Shell(int globalShell) const;

	/**
   * Method pour connaitre le nombre de contractions pour une fonction de base
	 */
	int getNbrOfContractions4BasisFunction(int globalBasisFunction) const;


	/**
	 * Method Get for the molecule that is attach.
	 */
	const molecule * getMolecule(void) const;

	/**
	 * Method GET for the total number of shells that contains the atom item.
	 */
	int getNbrOfBasisFunctions4Shell(int atom , int shell) const;

	/**
	 * Method GET for the total number of shells that contains the atom item.
	 */
	int getNbrOfBasisFunctions4Shell(int shell) const;
	
	/**
	 * Method GET for the total number of shells that contains the atom item.
	 */
	int getNbrOfShells4Atom(int atom) const;

	/**
	 * Method GET to find the basis function position from a global number of 
	 * basis function.
	 */
	basisFunctionPosition getShellPosition4Global(int shellNumber) const;
	
	/**
	 * Method GET for the shells numerotation table.
	 */
	const containor<int> & getShells4Atoms(void) const;

	/**
	 * Method GET for the number of atoms.
	 */
	int getTotalNbrOfAtoms(void) const;

	/**
	 * Method GET for the total number of basis functions.
	 */
	int getTotalNbrOfBasisFunctions(void) const;

	/**
	 * Method GET for the total number of shells.
	 */
	int getTotalNbrOfShells(void) const;

	/**
	 * Method GET for the contraction og a specific shell.
	 */
	const contractions & getContractions4Shell(int shell) const;

	/**
		* Method GET for the contraction og a specific shell.
	 */
	const contractions & getContractions4BasisFunction(int basisFunction) const;
	

	/**
	 * Method GET for the contraction og a specific shell.
	 */
	const shellType & getShellType4Shell(int shell) const;
	
	/**
	 * Method GET for the contraction og a specific shell.
	 */
	const shellType & getShellType4BasisFunction(int basisFunction) const;
	

	/**
	 * Method SET for the basis functions numerotation table.
	 */
	void setBasisFunctions4Shells(const containor<int> & basisFunctions4Shells);

	/**
	 * Method SET for the basis functions numerotation table.
	 */
	void setShells4Atoms(const containor<int> & shells4Atoms);

	/**
	 * Method SET for the number of atoms.
	 */
	void setTotalNbrOfAtoms(int nbrOfAtoms);
};

#endif

