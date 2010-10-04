/* 
 * The model library of A.S.P.I.C. 
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
#ifndef _RESTRICTED_HARTREE_FOCK_
#define _RESTRICTED_HARTREE_FOCK_

#include "molecularSystem.h"
#include <matrixSymetric.h>
#include <matrixFull.h>
/**
 * Generic class to manipulate the hartree fock model.
 */

// Modification F. Legoll Nov 06
//class restrictedHartreeFock : public molecularSystem
class restrictedHartreeFock : public virtual molecularSystem
{
private:

protected:

public:
	
	/**
	 * Constructor of the class.
	 */
	restrictedHartreeFock(void);
	
	/**
	 * Destructor of the class.
	 */
	virtual ~restrictedHartreeFock(void);
	
	/**
	 * This is a pure virtual method to access the density matrix.
	 */
	virtual const matrixSymetric & getDensity(void) const = 0;

	/**
	 * Virtual method to set the coefficients.
	 */
	virtual void setCoefficients(const matrixFull & coefficients) = 0;

	/**
	 * This method compute the Energy of the molecule.
	 */
	double getElectronsEnergy(const matrixSymetric & density) const;	
	
	/**
	 * This method compute the Energy of the molecule.
	 */
	double getElectronsEnergy(void) const;	

	/**
	 * This method compute the total energy of the molecule.
	 */
	double getEnergy(const matrixSymetric & density) const;	
	
	/**
	 * This method compute the total energy of the molecule.
	 */
	double getEnergy(void) const;	
	
	/**
	 * This method computes the fock matrix.
	 */
	matrixSymetric getFockMatrix(void) const;

	/**
	 * This method computes the fock matrix.
	 */
	matrixSymetric getFockMatrix(const matrixSymetric & denisty) const;

	/**
	 * This method computes the fock matrix.
	 */
	matrixSymetric applyTensor4RHF(const matrixSymetric & density) const;

};



#endif
