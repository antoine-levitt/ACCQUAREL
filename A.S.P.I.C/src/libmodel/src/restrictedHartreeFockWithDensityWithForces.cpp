/* 
* The model library of the A.S.P.I.C. 
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
#include "restrictedHartreeFockWithDensityWithForces.h"
//#include <generalizedEigenSolver.h>
//#include <vector.h>


#include <generalizedEigenSolver.h>
#include "restrictedHartreeFock.h"
#include <matrixSymetric.h>
#include <vector.h>

#include "restrictedHartreeFockWithDensity.h"


//////////////////////////////////////////////////////////////////////////////////////////
// constructor of the class.
//////////////////////////////////////////////////////////////////////////////////////////
restrictedHartreeFockWithDensityWithForces::restrictedHartreeFockWithDensityWithForces(void)
: restrictedHartreeFockWithDensity()
{
	;
}

//////////////////////////////////////////////////////////////////////////////////////////
// destructor of the class.
//////////////////////////////////////////////////////////////////////////////////////////
restrictedHartreeFockWithDensityWithForces::~restrictedHartreeFockWithDensityWithForces(void)
{
	clear();
}

//////////////////////////////////////////////////////////////////////////////////////////
// clear method of the class.
//////////////////////////////////////////////////////////////////////////////////////////
void restrictedHartreeFockWithDensityWithForces::clear(void)
{
	restrictedHartreeFockWithDensity::clear();
	molecularSystemWithForces::clearNewMembers();
	
}

//////////////////////////////////////////////////////////////////////////////////////////
// the method to compute the energy weighted density matrix.
//////////////////////////////////////////////////////////////////////////////////////////
const matrixSymetric restrictedHartreeFockWithDensityWithForces::ComputeDensityE(void) const
{

  int nbrOfBasisFunctions = getNbrOfBasisFunctions();
  int nbrOfElectronsdiv2 = getNbrOfElectrons() /2;
  matrixSymetric DE(nbrOfBasisFunctions);

  matrixFull coefficients = getCoefficients();
  containor<double> EnergyLevels = getEigenValues();
  
  for(int i=0 ; i < nbrOfBasisFunctions ; i++) {
	for(int j=i ; j < nbrOfBasisFunctions ; j++) {
	  DE(i,j) = 0.;
	  for(int k=0 ; k < nbrOfElectronsdiv2 ; k++) {
		// le 2 a cause du spin: cf la definition de la matrice de Fock: c'est 1/2 de la matrice reele
		// en tenant compte du spin
		DE(i,j) += 2.*coefficients(i,k) * coefficients(j,k) * EnergyLevels[k];
	  }
	} 
  }
  
  return DE;
}

////////////////////////////////////////////////////////////////////////////////////////////
// This is the method SET molecularSystem that is reimplemented from here
////////////////////////////////////////////////////////////////////////////////////////////
void restrictedHartreeFockWithDensityWithForces::setMolecularSystem(const molecule & mol , bool showProgress)
{

  restrictedHartreeFockWithDensity::setMolecularSystem(mol,showProgress);
  molecularSystemWithForces::setMolecularSystemWithForces(mol,showProgress);

}

