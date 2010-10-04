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
#ifndef _MOLECULAR_SYSTEM_
#define _MOLECULAR_SYSTEM_

#include <matrixSymetric.h>
#include <molecule.h>
#include <tensorSymetric.h>

/**
 * Class that build the basis matrixes for a given molecule.
 */
class molecularSystem
{

private:
	
	/**
	 * The coulomb tensor.
	 */
	tensorSymetric Coulomb;
	
	/**
	 * The time used to compute the coulomb tensor in seconds.
   */
	double CoulombComputationTime;

	/**
	 * The hamiltonian matrix.
	 */
	matrixSymetric Hamilton;
	
	/**
	 * The time used to compute the hamilton matrix in seconds.
   */
	double HamiltonComputationTime;

	/**
	 * The number of electrons that are in the molecular system.
	 */
	int NbrOfElectrons;

	/**
	 * The Enegrgy of the nuclei Energy interaction.
	 */
	double NucleiEnergy;
	
	/**
	 * The overlap matrix.
	 */
	matrixSymetric Overlap;
	
	/**
	 * The time used to compute the overlap matrix in seconds.
   */
	double OverlapComputationTime;
	
	

protected:
	
	/**
	 * Method that computes the coulomb tensor for a given molecule.
	 *
	 * @warning no space is allocated at this point. The method only feels
	 * the coulomb tensor, space must be allocated with the setNbrOfBasisFunctions
	 * method.
	 */
	void computeCoulomb(const molecule & mol , bool showProgress = false);
	
	/**
	 * Method that computes the hamiltonian matrix for a given molecule.
	 *
	 * @warning no space is allocated at this point. The method only feels
	 * the coulomb tensor, space must be allocated with the setNbrOfBasisFunctions
	 * method.
	 */
	void computeHamilton(const molecule & mol , bool showProgress = false);
	
	/**
	 * Method that computes the nuclei energy.
	 */
	void computeNucleiEnergy(const molecule & mol);
	
	/**
	 * Method that computes the overlap matrix for a given molecule.
	 *
	 * @warning no space is allocated at this point. The method only feels
	 * the coulomb tensor, space must be allocated with the setNbrOfBasisFunctions
	 * method.
	 */
	void computeOverlap(const molecule & mol , bool showProgress = false);
	
	/**
	 * Method to set the number of basis functions.
	 *
	 * This method allocates spaces for the matrixes and the tensor.
	 * 
	 * @param nbrOfBasisFunctions the number of basis functions that 
	 * the molecular system contains.
	 */
	void setNbrOfBasisFunctions(const int & nbrOfBasisFunctions);
	
	/**
	 * Method to set the number of electrons. 
	 */
	void setNbrOfElectrons(const int & nbrOfElectrons);
	
public:
	
	/**
	 * Default constrctor.
	 */
	molecularSystem(void);
	
	/**
	 * Destructor.
	 */
	virtual ~molecularSystem(void);
	
	/**
	 * Method that applies the coulomb tensor to a symetric matrix.
	 */
	matrixSymetric applyCoulomb(const matrixSymetric & density) const;
	
	/**
	 * Method that applies the exchange tensor to a symetric matrix.
	 */
	matrixSymetric applyExchange(const matrixSymetric & density) const;
	
	/**
	 * The clear method that frees all the space for the object.
	 */
	void clear(void);
	
	/**
	 * Method to know if the object is currently empty.
	 *
	 * @return true if empty false if not.
	 */
	bool empty(void) const;
	
	/**
	 * The GET method for the Coulomb tensor.
	 *
	 * @return a constant reference to the coulomb tensor.
	 */
	const tensorSymetric & getCoulomb(void) const;


	/**
	 * The GET method for the coulomb computation time.
	 *
	 * @return the coulomb computation time.
	 */
	const double & getCoulombComputationTime(void) const;

	/**
	 * The GET method for the Hamilton matrix.
	 *
	 * @return a constant reference to the coulomb tensor.
	 */
	const matrixSymetric & getHamilton(void) const;
	
	/**
	 * The GET method for the hamilton computation time.
	 *
	 * @return the hamilton computation time.
	 */
	const double & getHamiltonComputationTime(void) const;

	/**
	 * Method get for the number of basis functions.
	 *
	 * @return the number of basis functions for the given molecular system.
	 */
	const int & getNbrOfBasisFunctions(void) const;
	
	/**
	 * Method get for the number of electrons.
	 *
	 * @return the number of electron for the given molecular system.
	 */
	const int & getNbrOfElectrons(void) const;
	
	/**
	 * Method get for the nuclei energy
	 */
	const double & getNucleiEnergy(void) const;
	
	
	/**
	 * The GET method for the Overlap matrix.
	 *
	 * @return a constant reference to the overlap matrix.
	 */
	const matrixSymetric & getOverlap(void) const;
	
	/**
	 * The GET method for the overlap computation time.
	 *
	 * @return the overlap computation time.
	 */
	const double & getOverlapComputationTime(void) const;

	/**
	 * This method constructs the matrixes and the tensor 
	 * for a given molecule.
	 *
	 * @param the molecule for which the matrixes and tensor will be build.
	 */
	void setMolecularSystem(const molecule & mol, bool ShowProgress = false);
};

#endif
