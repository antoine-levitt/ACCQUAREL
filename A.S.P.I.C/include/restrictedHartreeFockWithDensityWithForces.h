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
#ifndef _RESTRICTED_HARTREE_FOCK_WITH_DENSITY_WITH_FORCES_
#define _RESTRICTED_HARTREE_FOCK_WITH_DENSITY_WITH_FORCES_

//#include <generalizedEigenSolver.h>
#include "restrictedHartreeFockWithDensity.h"
#include "restrictedHartreeFockWithForces.h"
#include <matrixSymetric.h>
//#include <vector.h>

/**
 * Class that implements the Hartree Fock model with the storage of the
 * density matrix and forces computations.
 */

class restrictedHartreeFockWithDensityWithForces : public restrictedHartreeFockWithForces, public restrictedHartreeFockWithDensity
{

private:

protected:

public:
	
	/**
	 * Constructor of the class restricted Hartree Fock model with density storage and forces computation.
	 */
	restrictedHartreeFockWithDensityWithForces(void);
	
	/**
	 * Destructor of the restricted Hartree Fock model with density storage and forces computation.
	 */
	virtual ~restrictedHartreeFockWithDensityWithForces(void);

	/**
	 * The clear method of the class.
	 */
	void clear(void);
	
	/**
	 * The method to compute the energy weighted density matrix.
	 */
	virtual const matrixSymetric ComputeDensityE(void) const;
	
	/**
	 * Method that set the molecular system for a restricted Hartree Fock
	 * Object. We mainly reuse the RHF_WithDensity class function
	 */
	void setMolecularSystem(const molecule & mol , bool showProgress = false);
	
};

#endif
