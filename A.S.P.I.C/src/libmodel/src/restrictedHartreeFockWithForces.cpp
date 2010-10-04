/* 
 * The model library of A.S.P.I.C. 
 * Class written by Frederic Legoll Nov 06
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
#include "restrictedHartreeFockWithForces.h"
#include <moleculeMapper.h>

//////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor of the class. 
//////////////////////////////////////////////////////////////////////////////////////////////////
restrictedHartreeFockWithForces::restrictedHartreeFockWithForces(void)
: molecularSystemWithForces()
{
	;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Destructor of the class. 
//////////////////////////////////////////////////////////////////////////////////////////////////
restrictedHartreeFockWithForces::~restrictedHartreeFockWithForces(void)
{
	;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// This method applies derivative of both coulomb tensors (the coulomb and the exchange) to a density matrix
// The spin is taken into account
//////////////////////////////////////////////////////////////////////////////////////////////////
matrixSymetric restrictedHartreeFockWithForces::applyDerTensor4RHF(const molecule & mol, const matrixSymetric & density, int atom, int dim)
{
  assert(density.getNbrOfRows() == getNbrOfBasisFunctions());
  
  matrixSymetric result;
  // check whether DerCoulomb is up to date
  int DerAtomCoulomb_,DerCoordCoulomb_;
  getDerInfo(DerAtomCoulomb_,DerCoordCoulomb_);
  if ((DerAtomCoulomb_ != atom) || (DerCoordCoulomb_ != dim)) {
	computeDerCoulomb(mol,atom,dim);
  }
  result = 2. * applyDerCoulomb(density);
  result -= applyDerExchange(density);
  return result;
	
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// This method computes the RHF Forces.
//////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull restrictedHartreeFockWithForces::getElectronsForces(const molecule & mol, const matrixSymetric & density)
{

  int nbrOfBasisFunctions = getNbrOfBasisFunctions();
  assert(density.getNbrOfRows() == nbrOfBasisFunctions);
  
  int NbrOfAtoms_ = getNbrOfAtoms();
  matrixFull Forces_(NbrOfAtoms_,3);
  Forces_.setCoefficients(0.);

  matrixSymetric energyMatrix;
  matrixSymetric correctOverlap;

  // In the Spinless HF model, dW_elec/dlambda = dh/dlambda::D + 1/2 D::dA/dlambda::D - Tr(D_E dS/dlambda)
  // with A = A^+ - A^-
  //
  // In the RHF, dW_elec/dlambda = 2 dh/dlambda::D + 2 D::dA^+/dlambda::D - D::dA^-/dlambda::D - Tr(D_E dS/dlambda)
  // The quantity 2 D::dA^+/dlambda - D::dA^-/dlambda has been computed with applyDerTensor4RHF
  
  // we need the matrix D_E (= C E C^*), which will be the energy weighted density matrix
  matrixSymetric densityE = ComputeDensityE();
  double value;
  atom::distanceUnit unit_;
  
  for (int atom=0; atom < NbrOfAtoms_ ; atom++) {
	for (int dim=0; dim<3; dim++) {
	  // we compute the force on the atom "atom" in the direction "dim"
	  
	  energyMatrix = 2*getDerHamilton(mol,atom,dim);
	  energyMatrix += applyDerTensor4RHF(mol,density,atom,dim);
	  correctOverlap = getDerOverlap(mol,atom,dim);
	  
	  value = 0.;
	  for(int i=0; i < nbrOfBasisFunctions ; i++) {
		for(int k=0; k < nbrOfBasisFunctions ; k++) {
		  value += energyMatrix(i,k) * density(k,i);
		  value -= correctOverlap(i,k) * densityE(k,i);
		}
	  }
	  // we now take into account the fact that the distance unit may not be the atomic unit:
	  unit_ = mol.getAtom(atom).getDistanceUnit();
	  if (unit_ != atom::ATOMIC_UNIT) {
		if(unit_ == atom::ANGSTROEM) {
		  value = value*(atom::AtomicUnit2Angstroem);
		} else {
		  assert(0);
		}
	  }
	  Forces_.setCoefficient(atom,dim,-value);
	}
  }
	
  return Forces_;
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// This method computes the RHF Forces
//////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull restrictedHartreeFockWithForces::getElectronsForces(const molecule & mol)
{
	return getElectronsForces(mol,getDensity());
}	

//////////////////////////////////////////////////////////////////////////////////////////////////
// This method computes the total forces
//////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull restrictedHartreeFockWithForces::getForces(const molecule & mol, const matrixSymetric & density)
{
  matrixFull A = getElectronsForces(mol,density);
  A.add(getNucleiForces());
  return A;
}	


//////////////////////////////////////////////////////////////////////////////////////////////////
// This method computes the total forces
//////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull restrictedHartreeFockWithForces::getForces(const molecule & mol)
{
	return getForces(mol,getDensity());
}	

