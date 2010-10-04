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
#ifndef _RESTRICTED_HARTREE_FOCK_WITH_FORCES_
#define _RESTRICTED_HARTREE_FOCK_WITH_FORCES_

#include "molecularSystemWithForces.h"
#include "restrictedHartreeFock.h"
#include <matrixSymetric.h>
#include <matrixFull.h>
/**
 * Generic class to manipulate the hartree fock model, with forces.
 */
class restrictedHartreeFockWithForces : public molecularSystemWithForces, public virtual restrictedHartreeFock
{
private:

protected:

public:
	
	/**
	 * Constructor of the class.
	 */
	restrictedHartreeFockWithForces(void);
	
	/**
	 * Destructor of the class.
	 */
	virtual ~restrictedHartreeFockWithForces(void);

	/**
	 * This method checks that derivatives of both coulomb tensors (coulomb and exchange) have been computed and then applies these derivatives to a density matrix
	 */
	matrixSymetric applyDerTensor4RHF(const molecule & mol, const matrixSymetric & density, int atom, int dim);
	
	/**
	 * This method computes the electronic forces on the molecule.
	 */
	matrixFull getElectronsForces(const molecule & mol, const matrixSymetric & density);

	/**
	 * This method computes the electronic forces on the molecule.
	 */
	matrixFull getElectronsForces(const molecule & mol);

	/**
	 * This method computes the total forces on the molecule
	 */
	matrixFull getForces(const molecule & mol, const matrixSymetric & density);
	
	/**
	 * This method computes the total forces on the molecule
	 */
	matrixFull getForces(const molecule & mol);
	
	/**
	 * This is a pure virtual method to compute the energy weighted density matrix
	 */
	virtual const matrixSymetric ComputeDensityE(void) const = 0;
	
};



#endif
