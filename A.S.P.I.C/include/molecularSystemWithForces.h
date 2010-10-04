/* 
 * The model library of the A.S.P.I.C. 
 * Class written by Frederic Legoll, Nov 06.
 *
 * Copyright (C) 2006 Frederic Legoll
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
#ifndef _MOLECULAR_SYSTEM_WITH_FORCES_
#define _MOLECULAR_SYSTEM_WITH_FORCES_

#include <matrixFull.h>
#include <molecule.h>
#include "molecularSystem.h"

/**
 * Class to add the forces to the class molecularSystem
 */
class molecularSystemWithForces: public virtual molecularSystem
{

private:
	
	/**
	 * The forces from nuclei interaction
	 */
     matrixFull NucleiForces;

	/**
	 * The derivative of the coulomb tensor wrt a nuclei coordinate.
	 */
	tensorSymetric DerCoulomb;
	
	/**
	 * The derivative of the hamiltonian matrix wrt a nuclei coordinate.
	 */
	matrixSymetric DerHamilton;
	
	/**
	 * The derivative of the overlap matrix wrt a nuclei coordinate.
	 */
	matrixSymetric DerOverlap;

	/**
	 * The number of atoms that are in the molecular system.
	 */
	int NbrOfAtoms;

	/**
	 * The number of the atom wrt which the derivate (of Hamilton, Coulomb, Overlap) has been computed
	 */
	int DerAtomHamilton, DerAtomOverlap, DerAtomCoulomb;

	/**
	 * The coordinate wrt which the derivate (of Hamilton, Coulomb, Overlap) has been computed
	 */
	int DerCoordHamilton, DerCoordOverlap, DerCoordCoulomb;
	
protected:

	/**
	 * Method to set the number of atoms. 
	 */
	void setNbrOfAtoms(const int & nbrOfAtoms);
	
	/**
	 * Method that computes the nuclei forces
	 */
	void computeNucleiForces(const molecule & mol);
	
	/**
	 * Method that computes the derivative of the overlap coeff (ij) wrt the nuclei p in the direction der_
	 */
	double computeOverlapDer(const molecule & mol, int i, int j, int p, ipoint<3> der_) const;
	
	/**
	 * Method to compute the derivative of the bielectronic integral (ijkl) wrt the nuclei coordinate (p,der_)
	 */
	double computeBielecDer(const molecule & mol, int i, int j, int k, int l, int p, ipoint<3> der_) const;
	
	/**
	 * Method that computes the derivative of the hamilton matrix wrt the nuclei coordinate (l,direction)
	 */
	void computeDerHamilton(const molecule & mol, int l, int direction);

	/**
	 * Method that computes the derivative of the overlap matrix wrt the nuclei coordinate (l,direction)
	 */
	void computeDerOverlap(const molecule & mol, int l, int direction);

	/**
	 * Method that computes the derivative of the coulomb tensor wrt the nuclei coordinate (l,direction)
	 */
	void computeDerCoulomb(const molecule & mol, int l, int direction);

	/**
	 * Method that applies the derivative of the coulomb tensor to a symetric matrix.
	 */
	matrixSymetric applyDerCoulomb(const matrixSymetric & density) const;
	
	/**
	 * Method that applies the derivative of the exchange tensor to a symetric matrix.
	 */
	matrixSymetric applyDerExchange(const matrixSymetric & density) const;
	
	
public:
	
	/**
	 * Default constructor.
	 */
	molecularSystemWithForces(void);
	
	/**
	 * Destructor.
	 */
	virtual ~molecularSystemWithForces(void);
	
	/**
	 * The clear method that frees all the space for the object.
	 */
	void clear(void);
	
	/**
	 * The clear method that frees the members that are new wrt to the molecularSystem class
	 */
	void clearNewMembers(void);

	/**
	 * Method GET for the number of atoms
	 */
	const int & getNbrOfAtoms(void) const;

	/**
	 * Method GET for the nuclei forces
	 */
	const matrixFull & getNucleiForces(void) const;

	/**
	 * The GET method for the derivative of Hamilton matrix.
	 */
	matrixSymetric & getDerHamilton(const molecule & mol, int atom, int dim);

	/**
	 * The GET method for the derivative of Overlap matrix.
	 */
	matrixSymetric & getDerOverlap(const molecule & mol, int atom, int dim);

	/**
	 * A method to know with respect to which atom/coord the derivative of Coulomb tensor has been lastly computed
	 */
	void getDerInfo(int & DerAtomCoulomb_, int & DerCoordCoulomb_) const;
	
	/**
	 * This method sets the correct size for the matrices and tensors 
	 */
	void setMolecularSystemWithForces(const molecule & mol, bool ShowProgress = false);
};

#endif
