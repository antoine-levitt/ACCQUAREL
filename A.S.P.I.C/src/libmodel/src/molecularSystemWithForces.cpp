/* 
 * The model library of the A.S.P.I.C. 
 * Class written by Frederic Legoll
 *
 * Copyright (C) 2006  Frederic Legoll
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
#include "molecularSystemWithForces.h"
#include <moleculeMapper.h>
#include <coulombIntegral.h>
#include <potentialIntegral.h>
#include <monoElectronicIntegral.h>
#include <time.h>

///////////////////////////////////////////////////////////////////////////////////
// Constructor of the class molecule.
///////////////////////////////////////////////////////////////////////////////////
molecularSystemWithForces::molecularSystemWithForces(void)
:molecularSystem() , 
NucleiForces()
{
	;
}

///////////////////////////////////////////////////////////////////////////////////
// Destructor of the class molecule.
///////////////////////////////////////////////////////////////////////////////////
molecularSystemWithForces::~molecularSystemWithForces(void)
{
	;
}


///////////////////////////////////////////////////////////////////////////////////
// Method clear for the molecular system class
///////////////////////////////////////////////////////////////////////////////////
void molecularSystemWithForces::clear(void)
{
  molecularSystem::clear();
  clearNewMembers();
}

///////////////////////////////////////////////////////////////////////////////////
// Method to clear the members that are new wrt to the MolecularSystem class
///////////////////////////////////////////////////////////////////////////////////
void molecularSystemWithForces::clearNewMembers(void)
{
  NucleiForces.clear();
  DerHamilton.clear();
  DerOverlap.clear();
  DerCoulomb.clear();
  
}

////////////////////////////////////////////////////////////////////////////////////
// Method that computes the nuclei forces
////////////////////////////////////////////////////////////////////////////////////
void molecularSystemWithForces::computeNucleiForces(const molecule & mol)
{

  int zk, zl;
  dpoint<3> xk , xl;
  double dist3,temp;
  	
  NucleiForces.setCoefficients(0.);
  int nbrOfAtoms_ = mol.getNbrOfAtoms();
  assert(nbrOfAtoms_ == NbrOfAtoms);

  for(int k=0 ; k < nbrOfAtoms_ ; k++) {
	xk = mol.getAtom(k).getPosition(atom::ATOMIC_UNIT);
	zk = mol.getAtom(k).getNbrOfProtons();	
	for(int l=k+1 ; l < nbrOfAtoms_ ; l++) {
	  xl = mol.getAtom(l).getPosition(atom::ATOMIC_UNIT);
	  zl = mol.getAtom(l).getNbrOfProtons();
	  dist3 = pow((xk-xl).norme_2(),3);
	  for (int dim=0; dim<3; dim++) {
		temp = zk*zl*(xk[dim] - xl[dim]) / dist3;
		NucleiForces.setCoefficient(k,dim,NucleiForces(k,dim)+temp);
		NucleiForces.setCoefficient(l,dim,NucleiForces(l,dim)-temp);
	  }
	}		
  }

  // we now take into account the fact that, for some atoms, the distance unit is not the atomic unit:
  atom::distanceUnit unit_;
  for(int k=0 ; k < nbrOfAtoms_ ; k++) {
	unit_ = mol.getAtom(k).getDistanceUnit();
	if (unit_ != atom::ATOMIC_UNIT) {
	  if(unit_ == atom::ANGSTROEM) {
		for (int dim=0; dim<3; dim++) {
		  NucleiForces.setCoefficient(k,dim,NucleiForces(k,dim)*(atom::AtomicUnit2Angstroem));
		}
	  } else {
		assert(0);
	  }
	}
  }
}


////////////////////////////////////////////////////////////////////////////////////
// Method to compute the derivative of an overlap wrt to a nuclei coord
////////////////////////////////////////////////////////////////////////////////////
double molecularSystemWithForces::computeOverlapDer(const molecule & mol, int i, int j, int p, ipoint<3> der_) const
{

  // value = d [\int Phi_i Phi_j] / dX_p^alpha
  // Notation: X_p^alpha is the coord #alpha of the nuclei #p 
  // the derivative direction is given by der_
  
  moleculeMapper map;
  monoElectronicIntegral overlap;
  gaussianBasisFunction phi_a , phi_b;

  map.attach(mol);

  ipoint<3> der0(0,0,0);
  		
  phi_a = map.getBasisFunction4Global(i);
  // the number of the atom to which the basis function i is attached
  int sigma_i = map.getBasisFunctionPosition4Global(i).Atom;
  phi_b = map.getBasisFunction4Global(j);
  // the number of the atom to which the basis function j is attached
  int sigma_j = map.getBasisFunctionPosition4Global(j).Atom;

  double value = 0.;
  if (p == sigma_i) {
	value -= overlap.getOverlapValue(phi_a,der_,phi_b,der0);
  }
  if (p == sigma_j) {
	value -= overlap.getOverlapValue(phi_a,der0,phi_b,der_);
  }

  return value;

}



////////////////////////////////////////////////////////////////////////////////////
// Method to compute the derivative of a bielectronic integral with respect to a nuclei coordinate
////////////////////////////////////////////////////////////////////////////////////
double molecularSystemWithForces::computeBielecDer(const molecule & mol, int i, int j, int k, int l, int p, ipoint<3> der_) const
{
  
  // valeur = d [\int Phi_i(x) Phi_j(x) Phi_k(x') Phi_l(x') / |x - x'| dx dx'] / dX_p^alpha
  // Notation: X_p^alpha is the coord #alpha of the nuclei #p 

  moleculeMapper map;
  coulombIntegral coulomb;
  gaussianBasisFunction phi_a , phi_b, phi_c, phi_d;

  map.attach(mol);

  ipoint<3> der0(0,0,0);
  		
  phi_a = map.getBasisFunction4Global(i);
  // the number of the atom to which the basis function i is attached
  int sigma_i = map.getBasisFunctionPosition4Global(i).Atom;
  phi_b = map.getBasisFunction4Global(j);
  // the number of the atom to which the basis function j is attached
  int sigma_j = map.getBasisFunctionPosition4Global(j).Atom;
  phi_c = map.getBasisFunction4Global(k);
  // the number of the atom to which the basis function k is attached
  int sigma_k = map.getBasisFunctionPosition4Global(k).Atom;
  phi_d = map.getBasisFunction4Global(l);
  // the number of the atom to which the basis function l is attached
  int sigma_l = map.getBasisFunctionPosition4Global(l).Atom;

  double value = 0.;
  if (p == sigma_i) {
	value -= coulomb.getCoulombValue(phi_a,der_,phi_b,der0,phi_c,der0,phi_d,der0);
  }
  if (p == sigma_j) {
	value -= coulomb.getCoulombValue(phi_a,der0,phi_b,der_,phi_c,der0,phi_d,der0);
  }
  if (p == sigma_k) {
	value -= coulomb.getCoulombValue(phi_a,der0,phi_b,der0,phi_c,der_,phi_d,der0);
  }
  if (p == sigma_l) {
	value -= coulomb.getCoulombValue(phi_a,der0,phi_b,der0,phi_c,der0,phi_d,der_);
  }

  return value;
  
}


////////////////////////////////////////////////////////////////////////////////////
// Method to compute the derivative of the Hamilton matrix wrt a nuclei coordinate
////////////////////////////////////////////////////////////////////////////////////
void molecularSystemWithForces::computeDerHamilton(const molecule & mol, int l, int direction)
{

  // we want to compute the derivative of h with respect to X_l^alpha
  // Notation: X_p^alpha is the coord #alpha of the nuclei #p
  // the direction alpha is encoded in der_
  
  moleculeMapper map;
  map.attach(mol);
  
  assert(l >= 0);
  assert(l < mol.getNbrOfAtoms());
  assert (direction >= 0);
  assert (direction <3);

  ipoint<3> der_;
  if (direction == 0) {
	der_.setDatas(1,0,0);
  }
  if (direction == 1) {
	der_.setDatas(0,1,0);
  }
  if (direction == 2) {
	der_.setDatas(0,0,1);
  }

  potentialIntegral potential;
  monoElectronicIntegral overlap;
  gaussianBasisFunction phi_a , phi_b;

  ipoint<3> der0(0,0,0);
  double dummy;

  int nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();

  int sigma_i, sigma_j;

  dpoint<3> X_bar_l = mol.getPosition4Atom(l);
  atom::distanceUnit unit_ = mol.getDistanceUnit4Atom(l);

  time_t start,end;
  time(&start);
  for (int i=0;i < nbrOfBasisFunctions ; i++) {	
    phi_a = map.getBasisFunction4Global(i);
	// the number of the atom to which the basis function i is attached
	sigma_i = map.getBasisFunctionPosition4Global(i).Atom;

	// DerHamilton is symetric: only consider j>=i:
	for(int j=i ; j < nbrOfBasisFunctions ; j++) {
	  phi_b = map.getBasisFunction4Global(j);
	  // the number of the atom to which the basis function j is attached
	  sigma_j = map.getBasisFunctionPosition4Global(j).Atom;

	  dummy = 0.;
	  
	  // taking into account kinetic energy
	  if (l == sigma_i) {
		dummy -= 0.5*(overlap.getKineticValue(phi_a,der_,phi_b,der0));
	  }
	  if (l == sigma_j) {
		dummy -= 0.5*(overlap.getKineticValue(phi_a,der0,phi_b,der_));
	  }

	  // we now take into account potential energy
	  dummy -= mol.getNbrOfProtons4Atom(l) * potential.getPotentialValue(phi_a,der_,phi_b,der0,X_bar_l,unit_);
	  dummy -= mol.getNbrOfProtons4Atom(l) * potential.getPotentialValue(phi_a,der0,phi_b,der_,X_bar_l,unit_);

	  // the last term coming from the potential energy
	  if (sigma_i == l) {
		for (int atom=0; atom < mol.getNbrOfAtoms(); atom++) {
		  dummy += mol.getNbrOfProtons4Atom(atom) * potential.getPotentialValue(phi_a,der_,phi_b,der0,mol.getPosition4Atom(atom),mol.getDistanceUnit4Atom(atom));
		}
		if (sigma_j == l) {
		  for (int atom=0; atom < mol.getNbrOfAtoms(); atom++) {
			dummy += mol.getNbrOfProtons4Atom(atom) * potential.getPotentialValue(phi_a,der0,phi_b,der_,mol.getPosition4Atom(atom),mol.getDistanceUnit4Atom(atom));
		  }
		}
	  }
	  if ((sigma_j == l) && (sigma_i != l)) {
		for (int atom=0; atom < mol.getNbrOfAtoms(); atom++) {
		  dummy += mol.getNbrOfProtons4Atom(atom) * potential.getPotentialValue(phi_a,der0,phi_b,der_,mol.getPosition4Atom(atom),mol.getDistanceUnit4Atom(atom));
		}
	  }
	  DerHamilton.setCoefficient(i,j,dummy);
	}
  }
  time(&end);
  //  cout<<"# Time to compute derivative of Hamiltonian matrix= "<<difftime(end,start)<<" seconds"<<endl;

  DerAtomHamilton = l;
  DerCoordHamilton = direction;

}	


////////////////////////////////////////////////////////////////////////////////////
// Method to compute the derivative of the Overlap matrix wrt a nuclei coordinate
////////////////////////////////////////////////////////////////////////////////////
void molecularSystemWithForces::computeDerOverlap(const molecule & mol, int l, int direction)
{

  // we want to compute the derivative of S with respect to X_l^alpha
  // Notation: X_p^alpha is the coord #alpha of the nuclei #p
  // the direction alpha is encoded in der_
  
  moleculeMapper map;
  map.attach(mol);
  
  assert(l >= 0);
  assert(l < mol.getNbrOfAtoms());
  assert (direction >= 0);
  assert (direction <3);

  ipoint<3> der_;
  if (direction == 0) {
	der_.setDatas(1,0,0);
  }
  if (direction == 1) {
	der_.setDatas(0,1,0);
  }
  if (direction == 2) {
	der_.setDatas(0,0,1);
  }
  
  int nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();

  double dummy;

  time_t start,end;
  time(&start);
  for (int i=0;i < nbrOfBasisFunctions ; i++) {	
	// DerOverlap is symetric: only consider j>=i:
	for(int j=i ; j < nbrOfBasisFunctions ; j++) {

	  dummy = computeOverlapDer(mol,i,j,l,der_);
	  DerOverlap.setCoefficient(i,j,dummy);
	}
  }
  time(&end);
  //  cout<<"# Time to compute derivative of Overlap matrix = "<<difftime(end,start)<<" seconds"<<endl;
  
  DerAtomOverlap = l;
  DerCoordOverlap = direction;
  
}	


////////////////////////////////////////////////////////////////////////////////////
// Method to compute the derivative of the Coulomb tensor wrt a nuclei coordinate
////////////////////////////////////////////////////////////////////////////////////
void molecularSystemWithForces::computeDerCoulomb(const molecule & mol, int p, int direction)
{

  // we want to compute the derivative of A = (ij|kl) with respect to X_p^alpha
  // Notation: X_p^alpha is the coord #alpha of the nuclei #p
  // the direction alpha is encoded in der_
  
  moleculeMapper map;
  map.attach(mol);
  
  assert(p >= 0);
  assert(p < mol.getNbrOfAtoms());
  assert (direction >= 0);
  assert (direction <3);

  int nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();

  double dummy;

  ipoint<3> der_;
  if (direction == 0) {
	der_.setDatas(1,0,0);
  }
  if (direction == 1) {
	der_.setDatas(0,1,0);
  }
  if (direction == 2) {
	der_.setDatas(0,0,1);
  }
  
	
  // on the loops (i,j,k,l), we take into considerations the symetry of DerCoulomb, which are exactly the same as the symetries of Coulomb
  time_t start,end;
  time(&start);
  for (int i=0; i < nbrOfBasisFunctions; i++) {	
	for(int j=i; j < nbrOfBasisFunctions; j++) {
	  for(int k=i; k < nbrOfBasisFunctions; k++) {
		for(int l=k; l < nbrOfBasisFunctions; l++) {
		  if (k==i && l<j) {
			continue;
		  } 

		  dummy = computeBielecDer(mol,i,j,k,l,p,der_);
		  DerCoulomb.setCoefficient(i,j,k,l,dummy);
		} // fin du for l.
	  } // fin du for k.	
	} // fin du for j.
  } // fin du for i.
  time(&end);
  //  cout<<"# Time to compute derivative of Coulomb matrix = "<<difftime(end,start)<<" seconds"<<endl;

  DerAtomCoulomb = p;
  DerCoordCoulomb = direction;
  
}

////////////////////////////////////////////////////////////////////////////////////
// Method that applies the derivative of coulomb tensor to a density matrix.
////////////////////////////////////////////////////////////////////////////////////
matrixSymetric molecularSystemWithForces::applyDerCoulomb(const matrixSymetric & density) const
{

  assert(density.getNbrOfRows() == getNbrOfBasisFunctions());
  
  matrixSymetric DerCoulombMatrix;
	
  if (DerCoulomb.empty()) {
	return DerCoulombMatrix;
  } 
	
  int nbrOfRows = getNbrOfBasisFunctions();
  double value;
	
  DerCoulombMatrix.setMatrixSize(getNbrOfBasisFunctions());
	
  for(int i=0 ; i < nbrOfRows ; i++) {
	for(int j=i ; j < nbrOfRows ; j++) {
	  value = 0.;
			
	  for(int k=0 ; k < nbrOfRows ; k++) {
		value += DerCoulomb.getCoefficient(i,j,k,k) * density(k,k);
		for(int l=k+1 ; l < nbrOfRows ; l++) {
		  value += 2*DerCoulomb.getCoefficient(i,j,k,l) * density.getCoefficient(k,l);
		}		
	  }
	  
	  DerCoulombMatrix.setCoefficient(i,j,value);
	}		
  }
  
  return DerCoulombMatrix;
}

////////////////////////////////////////////////////////////////////////////////////
// Method that applies the derivative of exchange tensor to a density matrix.
////////////////////////////////////////////////////////////////////////////////////
matrixSymetric molecularSystemWithForces::applyDerExchange(const matrixSymetric & density) const
{
	
  assert(density.getNbrOfRows() == getNbrOfBasisFunctions());
	
  matrixSymetric DerExchangeMatrix;
	
  if (DerCoulomb.empty()) {
	return DerExchangeMatrix;
  }
	
  int nbrOfRows = getNbrOfBasisFunctions();
  double value;
	
  DerExchangeMatrix.setMatrixSize(nbrOfRows);
	
  for(int i=0 ; i < nbrOfRows ; i++) {
	for(int j=i ; j < nbrOfRows ; j++) {
			
	  value = 0.;
	  for(int k=0 ; k < nbrOfRows ; k++) {
		for(int l=0 ; l < nbrOfRows ; l++) {
		  value += DerCoulomb.getCoefficient(i,k,j,l) * density.getCoefficient(k,l);
		}		
	  }
			
	  DerExchangeMatrix.setCoefficient(i,j,value);
	}		
  }
	
  return DerExchangeMatrix;
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the nuclei forces
///////////////////////////////////////////////////////////////////////////////////////
const matrixFull & molecularSystemWithForces::getNucleiForces(void) const
{
	return NucleiForces;
}

////////////////////////////////////////////////////////////////////////////////////////
// Method SET from a molecule.
////////////////////////////////////////////////////////////////////////////////////////
void molecularSystemWithForces::setMolecularSystemWithForces(const molecule & mol , bool showProgress)
{

  setMolecularSystem(mol,showProgress);
  
  int nbrOfBasisFunctions = mol.getNbrOfBasisFunctions();
  DerCoulomb.setTensorSize(nbrOfBasisFunctions);
  DerHamilton.setMatrixSize(nbrOfBasisFunctions);
  DerOverlap.setMatrixSize(nbrOfBasisFunctions);
  setNbrOfAtoms(mol.getNbrOfAtoms());
  // we do not compute the derivative of Coulomb, Hamilton, Overlap since there is no particular direction/nuclei
  DerAtomHamilton = -1;
  DerAtomOverlap = -1;
  DerAtomCoulomb = -1;
  DerCoordHamilton = -1;
  DerCoordOverlap = -1;
  DerCoordCoulomb = -1;

  NucleiForces.setMatrixSize(mol.getNbrOfAtoms(),3);
  computeNucleiForces(mol);  
  
}

////////////////////////////////////////////////////////////////////////////////////////
// Method SET for the number of atoms
////////////////////////////////////////////////////////////////////////////////////////
void molecularSystemWithForces::setNbrOfAtoms(const int & nbrOfAtoms)
{
	assert(nbrOfAtoms >= 0);
	
	NbrOfAtoms = nbrOfAtoms;
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the number of atoms
///////////////////////////////////////////////////////////////////////////////////////
const int & molecularSystemWithForces::getNbrOfAtoms(void) const
{
	return NbrOfAtoms;
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the derivative of Hamiltonian matrix.
///////////////////////////////////////////////////////////////////////////////////////
matrixSymetric & molecularSystemWithForces::getDerHamilton(const molecule & mol, int atom, int dim)
{
  // check that DerHamilton is up to date
  if ((atom != DerAtomHamilton) || (dim != DerCoordHamilton)) {
	computeDerHamilton(mol,atom,dim);
  }
  
  return DerHamilton;
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the derivative of Overlap matrix.
///////////////////////////////////////////////////////////////////////////////////////
matrixSymetric & molecularSystemWithForces::getDerOverlap(const molecule & mol, int atom, int dim)
{
  // check that DerOverlap is up to date
  if ((atom != DerAtomOverlap) || (dim != DerCoordOverlap)) {
	computeDerOverlap(mol,atom,dim);
  }
  
  return DerOverlap;
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the derivative information
///////////////////////////////////////////////////////////////////////////////////////
void molecularSystemWithForces::getDerInfo(int & DerAtomCoulomb_, int & DerCoordCoulomb_) const
{
  
  DerAtomCoulomb_   = DerAtomCoulomb;
  DerCoordCoulomb_  = DerCoordCoulomb;
}

