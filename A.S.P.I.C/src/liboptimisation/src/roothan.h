/* 
 * The optimisation library of the A.S.P.I.C. 
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

#ifndef _ROOTHAN_
#define _ROOTHAN_

#include <containor.h>
#include <restrictedHartreeFock.h>
#include <vector.h>

/**
 * This class performs a Roothan optimisation
 * on a restricted Hartree-Fock Model.
 */
class roothan4RestrictedHartreeFock
{
private:
	/**
	 * The number of iterations.
	 */
	int NbrOfIterations;

	/**
	 * The maximal number of iterations that shall be performed.
	 */
	int NbrOfIterationsMax;

	/**
	 * Those valriable stores the density matrix needed for 
	 * the roothan algorithm computation.
	 */
	matrixSymetric Densitys[2];

	/**
	 * The eigenvalues.
	 */
	containor<double> EigenValues[2];

	/**
	 * This contains the max value of the error to ensure convergence.
	 */
	double Max4Error;

protected:

	/**
	 * Copy constructor of the class.
	 *
	 * This constructor shall not be used. This implementation
	 * exists only to ensure he'll be never called.
	 */
	roothan4RestrictedHartreeFock(const roothan4RestrictedHartreeFock & roothan4rhf);

	/**
	 * Method that is used to perform the convergence test.
	 */
	bool convergenceHasBeenReached(void) const;

	/**
	 * Method to access one of the two density matrixes
	 * stored in the roothan object. This is the constant version.
	 *
	 * @param item if the item is odd then we return the first density
	 * else the second.
	 */
	const matrixSymetric & getDensity(const int & item) const;

	/**
	 * Method to access one of the eigenvalues
	 * stored in the roothan object. This is the constant version.
	 *
	 * @param item if the item is odd then we return the first density
	 * else the second.
	 */
	const containor<double> & getEigenValues(const int & item) const;

	/**
	 * Method to set one of the two density matrixes
	 * stored in the roothan object.
	 *
	 * @param item if the item is odd then we return the first density
	 * else the second.
	 *
	 * @param density the density matrix to set.
	 */
	 void setDensity(const int & item , const matrixSymetric & density);

	/**
	 * Method to set one of the two density matrixes
	 * stored in the roothan object.
	 *
	 * @param item if the item is odd then we return the first density
	 * else the second.
	 *
	 * @param density the density matrix to set.
	 */
	 void setDensity(const int & item , const containor<vector> & eigenVectors , const int & nbrOfElectrons);

	/**
	 * Method to set one of the eigenvalues
	 * stored in the roothan object. This is the constant version.
	 *
	 * @param item if the item is odd then we return the first density
	 * else the second.
	 */
	void setEigenValues(const int & item , const containor<double> & eigenValues);

	 /**
	 * Method to set the size for the denity matrix.
	 *
	 * @param nbrOfRows the number of row of the density matrix.
	 */
	void setNbrOfRows4Density(const int & nbrOfRows);

public:

	/**
	 * Default constructor of the class.
	 */
	roothan4RestrictedHartreeFock(void);

	/**
	 * Destructor of the class.
	 */
	~roothan4RestrictedHartreeFock(void);

	/**
	 * This method clears the current object.
	 */
	void clear(void);

	/**
	 * Method to know if we use an empty object.
	 *
	 * @return true if the current object is empty, false
	 * if not.
	 */
	bool empty(void) const;

	/**
	 * This is the acces method for the current matrix on which the 
	 * roothan algorithm has found.
   */
	const matrixSymetric & getDensity(void) const;

	/**
	 * This is the acces method for the current eigenvalues on which the 
	 * roothan algorithm has found.
   */
	const containor<double> & getEigenValues(void) const;

	/**
	 * Method that computes the Ground state of an
	 * restricted Hartree Fock Model using the roothan
	 * algorithm.
	 *
	 * The followin method changes the state given as an 
	 * argument in order to minimize the RHF energy of the given
	 * molecule.
	 */
	const matrixSymetric & getGroundStateDensity(const restrictedHartreeFock & rhfState , const bool & log = false);

 /**
  * Method to acces the maximal error value under which
	* the algorithm will be considered as convergent.
	*/
	const double & getMax4Error(void) const;

	/**
	 * This method is used to acces the number of iterations 
	 * that have been performed by the roothan algorithm.
   */
	const int & getNbrOfIterations(void) const;

	/**
	 * Method to get the size for the denity matrix.
	 *
	 * @return the number of row of the optimal density matrix.
	 */
	const int & getNbrOfRows4Density(void) const;

	/**
	 * This method is used to acces the maximal number of iterations 
	 * that shall be performed during the roothan loop.
   */
	const int & getNbrOfIterationsMax(void) const;

	/**
	 * The set method for the maximal value that is allowed on 
	 * the error.
	 */
	void setMax4Error(const double & max4Error);

	/**
	 * This method is used to set the maximal number of iterations
	 * that shall be performed during the roothan loop.
	 */
	void setNbrOfIterationsMax(const int & nbrOfIterationsMax);
};

#endif