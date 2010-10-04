/* 
 * The model library of the A.S.P.I.C. 
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
#ifndef _RESTRICTED_HARTREE_FOCK_WITH_DENSITY_
#define _RESTRICTED_HARTREE_FOCK_WITH_DENSITY_

#include <generalizedEigenSolver.h>
#include "restrictedHartreeFock.h"
#include <matrixSymetric.h>
#include <vector.h>

/**
 * Class that implements the Hartree Fock model with the storage of the
 * density matrix.
 */

// Modification F. Legoll Nov 06
//class restrictedHartreeFockWithDensity : public restrictedHartreeFock
class restrictedHartreeFockWithDensity : public virtual restrictedHartreeFock
{

private:

	/**
	 * The generalized eigen solver that converts the density matrix
	 * to the coefficients matrix.
	 */
	generalizedEigenSolver Density2Coefficients;
	
	/**
	 * The boolean tag that marks if the Density to Coefficients
	 * converter is up to date.
	 */
	bool Density2CoefficientsIsUp2Date;
	
	/**
	 * the density matrix.
	 */
	matrixSymetric Density;

protected:

	/**
   * Method to know if the density to coefficient
   * converter is up to date. 
   *
   * @return true if all is up to date false if the 
   * updateCoefficient method needs to be called.
   */
		bool density2CoefficientsIsUp2Date(void) const;
	
	/**
	 * Method that updates the Coefficients that
	 * are stored.
	 */
	void updateCoefficients(void);

public:
	
	/**
	 * Constructor of the class restricted Hartree Fock model with density storage.
	 */
	restrictedHartreeFockWithDensity(void);
	
	/**
	 * Destructor of the restricted Hartree Fock model wit density storage.
	 */
	virtual ~restrictedHartreeFockWithDensity(void);
	
	/**
	 * The clear method of the class.
	 */
	void clear(void);
	
	/**
		* The method that gives the coefficients for the current
	 * rhf state.
	 */
	matrixFull getCoefficients(void) const;

	/**
	 * The access method to the density matrix.
	 */
	virtual const matrixSymetric & getDensity(void) const;
	
	/**
		* The method that returns the Eigenvalues for the RHF 
	 * state.
	 */
	const containor<double> & getEigenValues(void) const;
	
	/**
	 * The method that returns the Eigenvalues for the RHF 
	 * state.
	 */
	const double & getEigenValue(const int & item) const;
	
	/**
		* The method that returns the Eigenvalues for the RHF 
	 * state.
	 */
	const vector & getEigenVector(const int & item) const;

	/**
	 * The method that returns the Eigenvalues for the RHF 
	 * state.
	 */
	const containor<vector> & getEigenVectors(void) const;
	
	/**
	 * Method to set the coefficient Matrix.
	 */
	virtual void setCoefficients(const matrixFull & coefficients);
	
	/**
	 * Method to set the density Matrix.
	 */
	virtual void setDensity(const matrixSymetric & density);
	
	/**
	 * Method that set the molecular system for a restricted Hartree Fock
	 * Object.
	 *
	 * This method inits the molecular system but in a mich more specific
	 * way to the restricted hartree fock model : This method verifys that 
	 * the number of electrons in the molecule is odd (needed to perform
	 * an RHF computation) and initialize the density matrix so that the object
	 * can be directly usefull from this point.
	 *
	 * @param mol the molecule for which the model shall be build.
	 */
	void setMolecularSystem(const molecule & mol , bool showProgress = false);
	
};

#endif
