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
#include "moleculeMapper.h"

/////////////////////////////////////////////////////////////////////////////////////////
// The default constructor.
/////////////////////////////////////////////////////////////////////////////////////////
moleculeMapper::moleculeMapper(void)
{
	clear();
}

/////////////////////////////////////////////////////////////////////////////////////////
// The destructor.
/////////////////////////////////////////////////////////////////////////////////////////
moleculeMapper::~moleculeMapper(void)
{
	clear();
}

/////////////////////////////////////////////////////////////////////////////////////////
// The Method to attach a molecule to the numérotation.
/////////////////////////////////////////////////////////////////////////////////////////
void moleculeMapper::attach(const molecule & mol) 
{

	// first step is to clean the object.
	clear();

	// we store a link to the molecule.
	setMolecule(&mol);

	// we compute the tables.
	computeTables();
}

/////////////////////////////////////////////////////////////////////////////////////////
// The Clear() method.
/////////////////////////////////////////////////////////////////////////////////////////
void moleculeMapper::clear(void)
{
	Molecule = NULL;
	BasisFunctions4Shells.clear();
	NbrOfAtoms = 0;
	Shells4Atoms.clear();
}

/////////////////////////////////////////////////////////////////////////////////////////
// The Method to compute the numerotation tables.
/////////////////////////////////////////////////////////////////////////////////////////
void moleculeMapper::computeBasisFunctionsTable(void)
{
	int i , j , nbrOfAtoms , nbrOfShells , nbrOfBasisFunctions , shell;

	nbrOfAtoms = getMolecule()->getNbrOfAtoms();
	nbrOfBasisFunctions = 0;
	shell = 0;

	nbrOfShells = getTotalNbrOfShells();

	if(nbrOfShells == 0) {
		return;
	}

	BasisFunctions4Shells.setSizes(nbrOfShells);
	
	for(i=0 ; i < nbrOfAtoms ; i++) {
		nbrOfShells = getMolecule()->getNbrOfShells4Atom(i);
		for(j=0 ; j < nbrOfShells ; j++) {
			nbrOfBasisFunctions += getMolecule()->getBasis4Atom(i).getNbrOfBasisFunctions4Shell(j);
			BasisFunctions4Shells[shell] = nbrOfBasisFunctions;
			shell++;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// The Method to compute the numerotation tables.
/////////////////////////////////////////////////////////////////////////////////////////
void moleculeMapper::computeShellsTable(void)
{
	int  i , nbrOfAtoms , nbrOfShells;

	nbrOfShells = 0;
	nbrOfAtoms = getTotalNbrOfAtoms();

	if(nbrOfAtoms == 0) {
		return;
	}

	// We compute the shells 4 atoms table.
	Shells4Atoms.setSizes(nbrOfAtoms);

	for(i=0 ; i < nbrOfAtoms ; i++) {
		nbrOfShells += getMolecule()->getNbrOfShells4Atom(i);
		Shells4Atoms[i] = nbrOfShells;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// The Method to compute the numerotation tables.
/////////////////////////////////////////////////////////////////////////////////////////
void moleculeMapper::computeTables(void)
{
	if(getMolecule() == NULL || getMolecule()->empty()) {
		return;
	}
	
	// we set the number of atoms.
	setTotalNbrOfAtoms(getMolecule()->getNbrOfAtoms());

	// we compute the shell table.
	computeShellsTable();
	
	// we compute the basis functions table.
	computeBasisFunctionsTable();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The copy method.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void moleculeMapper::copy(const moleculeMapper & map)
{
	if(map.empty()){
		clear();
	}

	setTotalNbrOfAtoms(map.getTotalNbrOfAtoms());
	setShells4Atoms(map.getShells4Atoms());
	setBasisFunctions4Shells(map.getBasisFunctions4Shells());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method to know if the object is empty.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool moleculeMapper::empty(void) const
{
	if(getTotalNbrOfAtoms() == 0 || getMolecule() == NULL) {
		return true;
	} else {
		return false;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for a basis function from a global numbering.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
gaussianBasisFunction moleculeMapper::getBasisFunction4Global(int basisFunctionNumber) const
{
	assert(basisFunctionNumber >=0);
	assert(basisFunctionNumber < getTotalNbrOfBasisFunctions());
		
	//gaussianBasisFunction function;
	basisFunctionPosition position;
		
	position = getBasisFunctionPosition4Global(basisFunctionNumber);


	return getMolecule()->getGaussianBasisFunction(position.Atom , position.Shell , position.BasisFunction);

	//function.setCenter(getMolecule()->getAtom(position.Atom).getPosition());
	//function.setContractions(getMolecule()->getBasis4Atom(position.Atom).getContractions4Shell(position.Shell));
	//function.setMonomeDegree(getMolecule()->getBasis4Atom(position.Atom).getMonomeDegree(position.Shell,position.BasisFunction));
	//return function;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for a basis function from a global numbering.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const dpoint<3> & moleculeMapper::getBasisFunctionCenter(int basisFunctionNumber) const
{
	assert(basisFunctionNumber >=0);
	assert(basisFunctionNumber < getTotalNbrOfBasisFunctions());
		
	basisFunctionPosition position;
		
	position = getBasisFunctionPosition4Global(basisFunctionNumber);

	return getMolecule()->getAtom(position.Atom).getPosition();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for a basis function from a global numbering.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const dpoint<3> & moleculeMapper::getShellCenter(int shell) const
{
	assert(shell >=0);
	assert(shell < getTotalNbrOfShells());
		
	basisFunctionPosition position;
		
	position = getShellPosition4Global(shell);

	return getMolecule()->getAtom(position.Atom).getPosition();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for a basis function from a global numbering.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const contractions &  moleculeMapper::getContractions4Shell(int shell) const
{
	assert(shell >=0);
	assert(shell < getTotalNbrOfShells());
	basisFunctionPosition position;
	position = getShellPosition4Global(shell);

	return getMolecule()->getBasis4Atom(position.Atom).getContractions4Shell(position.Shell);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for a basis function from a global numbering.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const contractions &  moleculeMapper::getContractions4BasisFunction(int basisFunction) const
{
	assert(basisFunction >=0);
	assert(basisFunction < getTotalNbrOfBasisFunctions());
	
	basisFunctionPosition position;
	position = getBasisFunctionPosition4Global(basisFunction);
	
	return getMolecule()->getBasis4Atom(position.Atom).getContractions4Shell(position.Shell);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for a basis function from a global numbering.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const shellType &  moleculeMapper::getShellType4BasisFunction(int basisFunction) const
{
	assert(basisFunction >=0);
	assert(basisFunction < getTotalNbrOfBasisFunctions());
	
	basisFunctionPosition position;
	position = getBasisFunctionPosition4Global(basisFunction);
	
	return getMolecule()->getBasis4Atom(position.Atom).getShellType(position.Shell);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for a basis function from a global numbering.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const shellType &  moleculeMapper::getShellType4Shell(int shell) const
{
	assert(shell >=0);
	assert(shell < getTotalNbrOfShells());
	
	basisFunctionPosition position;
	position = getShellPosition4Global(shell);
	
	return getMolecule()->getBasis4Atom(position.Atom).getShellType(position.Shell);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
basisFunctionPosition moleculeMapper::getBasisFunctionPosition4Global(int basisFunction) const
{
	basisFunctionPosition position;
	
	int shell;
	
	for(shell=0 ; shell < getTotalNbrOfShells() ; shell++) {
		if(basisFunction < BasisFunctions4Shells[shell]) {
			break;
		}
	}

	if(shell > 0) {
		basisFunction -= BasisFunctions4Shells[shell-1];
	}

	position = getShellPosition4Global(shell);
	position.BasisFunction = basisFunction;

	return position;
}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for the basis functions numerotation.
/////////////////////////////////////////////////////////////////////////////////////////
const containor<int> & moleculeMapper::getBasisFunctions4Shells(void) const
{
	return BasisFunctions4Shells;
}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for global number of a basis function.
/////////////////////////////////////////////////////////////////////////////////////////
int moleculeMapper::getGlobalBasisFunctionNumber(int atom , int shell , int basisFunction) const
{
	assert(atom >= 0);
	assert(atom < getTotalNbrOfAtoms());

	assert(shell >= 0);
	assert(shell < getNbrOfShells4Atom(atom));

	assert(basisFunction >=0);
	assert(basisFunction < getNbrOfBasisFunctions4Shell(atom,shell));

	if(atom > 0) {
		shell += Shells4Atoms[atom-1];
	} 
	
	if(shell > 0) {
		basisFunction += BasisFunctions4Shells[shell-1];
	}

	return basisFunction;
}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for global number of a shell.
/////////////////////////////////////////////////////////////////////////////////////////
int moleculeMapper::getGlobalShellNumber(int atom , int shell) const
{
	assert(atom >= 0);
	assert(atom < getTotalNbrOfAtoms());

	assert(shell >= 0);
	assert(shell < getNbrOfShells4Atom(atom));

	if(atom > 0) {
		shell += Shells4Atoms[atom-1];
	} 
		
	return  shell;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET for the molecule.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
const molecule * moleculeMapper::getMolecule(void) const
{
	return Molecule;
}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for number of basis functions in a shell
/////////////////////////////////////////////////////////////////////////////////////////
int moleculeMapper::getNbrOfBasisFunctions4Shell(int shell) const
{
	assert(shell >= 0);
	assert(shell < getTotalNbrOfShells());

	if(shell > 0) {
		return BasisFunctions4Shells[shell] - BasisFunctions4Shells[shell -1];	
	} else {
		return BasisFunctions4Shells[shell];	
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for number of basis functions in a shell
/////////////////////////////////////////////////////////////////////////////////////////
int moleculeMapper::getNbrOfBasisFunctions4Shell(int atom , int shell) const
{
	assert(atom >= 0);
	assert(atom < getTotalNbrOfAtoms());

	assert(shell >= 0);
	assert(shell < getNbrOfShells4Atom(atom));

	shell = getGlobalShellNumber(atom,shell);
	return getNbrOfBasisFunctions4Shell(shell);

}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for number of contraction for a given shell ( global numbering)
/////////////////////////////////////////////////////////////////////////////////////////
int moleculeMapper::getNbrOfContractions4Shell(int shell) const
{
	assert(shell >= 0);
	assert(shell < getTotalNbrOfShells());
	
	basisFunctionPosition position;
		
	position = getShellPosition4Global(shell);

	return getMolecule()->getNbrOfContractions4Shell(position.Atom , position.Shell);	
}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for number of contraction for a given shell ( global numbering)
/////////////////////////////////////////////////////////////////////////////////////////
int moleculeMapper::getNbrOfContractions4BasisFunction(int basisFunction) const
{
	assert(basisFunction >= 0);
	assert(basisFunction < getTotalNbrOfBasisFunctions());
	
	basisFunctionPosition position;
		
	position = getBasisFunctionPosition4Global(basisFunction);

	return getMolecule()->getNbrOfContractions4Shell(position.Atom , position.Shell);	
}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for number of shells in an atom
/////////////////////////////////////////////////////////////////////////////////////////
int moleculeMapper::getNbrOfShells4Atom(int atom) const
{
	assert(atom >= 0);
	assert(atom < getTotalNbrOfAtoms());

	if(atom == 0) {
		return Shells4Atoms[0];	
	} else {
		return Shells4Atoms[atom] - Shells4Atoms[atom -1];	
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method GET that gives the position of the local position of the shell and the atom for a global shell number.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
basisFunctionPosition moleculeMapper::getShellPosition4Global(int shell) const
{
	basisFunctionPosition position;
	
	int atom;
	
	for(atom = 0 ; atom < getTotalNbrOfAtoms() ; atom++) {
		if(shell < Shells4Atoms[atom]) {
			break;
		}
	}

	if(atom > 0) {
		shell -= Shells4Atoms[atom-1];
	}

	position.Atom = atom;
	position.BasisFunction = -1;
	position.Shell = shell;

	return position;
}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for the shell table numerotation.
/////////////////////////////////////////////////////////////////////////////////////////
const containor<int> & moleculeMapper::getShells4Atoms(void) const
{
	return Shells4Atoms;
}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for the number of atoms.
/////////////////////////////////////////////////////////////////////////////////////////
int moleculeMapper::getTotalNbrOfAtoms(void) const
{
	return NbrOfAtoms;
}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for the number of shells.
/////////////////////////////////////////////////////////////////////////////////////////
int moleculeMapper::getTotalNbrOfShells(void) const
{
	int nbrOfAtoms = getTotalNbrOfAtoms();
	
	if(nbrOfAtoms == 0) 
		return 0;

	return Shells4Atoms[nbrOfAtoms - 1];
}

/////////////////////////////////////////////////////////////////////////////////////////
// The GET method for the number of basis functions.
/////////////////////////////////////////////////////////////////////////////////////////
int moleculeMapper::getTotalNbrOfBasisFunctions(void) const
{
	int nbrOfShells = getTotalNbrOfShells();
	
	if(nbrOfShells == 0) 
		return 0;

	return BasisFunctions4Shells[nbrOfShells - 1];
}

/////////////////////////////////////////////////////////////////////////////////////////
// The SET method for the basis functions numerotation.
/////////////////////////////////////////////////////////////////////////////////////////
void moleculeMapper::setBasisFunctions4Shells(const containor<int> & basisFunctions4Shells)
{
	assert(getTotalNbrOfShells() == basisFunctions4Shells.getSizes());
	BasisFunctions4Shells = basisFunctions4Shells;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Method SET for the molecule.
/////////////////////////////////////////////////////////////////////////////////////////
void moleculeMapper::setMolecule(const molecule * mol)
{
	Molecule = mol;
}

/////////////////////////////////////////////////////////////////////////////////////////
// The SET method for the shells numerotation.
/////////////////////////////////////////////////////////////////////////////////////////
void moleculeMapper::setShells4Atoms(const containor<int> & shells4Atoms)
{
	assert(getTotalNbrOfAtoms() == shells4Atoms.getSizes());
	Shells4Atoms = shells4Atoms;
}

/////////////////////////////////////////////////////////////////////////////////////////
// The SET method for the number of atoms.
/////////////////////////////////////////////////////////////////////////////////////////
void moleculeMapper::setTotalNbrOfAtoms(int nbrOfAtoms)
{
	assert(nbrOfAtoms > 0);
	NbrOfAtoms = nbrOfAtoms;
}

