/* 
 * The martix library of the A.S.P.I.C. 
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

#ifndef _GENERALIZED_EIGEN_SOLVER_
#define _GENERALIZED_EIGEN_SOLVER_

#include <containor.h>
#include "vector.h"
#include "matrixFull.h"
#include "matrixSymetric.h"

/**
 * Classe pour résoudre un problème aux valeurs propres généralisées :
 * AX = lambda BX.
 *
 * A est une matrice ....
 * B est une matrice ....
 * 
 * Cette classe calcule (et stocke) toutes les valeurs propres (lambda) et
 * vecteurs propres (X) pour le problème.
 */
class generalizedEigenSolver
{
private:

	/**
	 * This matrix contains the eigen vectors for the problem.
	 */
	containor<vector> EigenVectors;

	/**
	 * This contains the eigne values.
	 */
	containor<double> EigenValues;

protected:

	/**
	 * This method allocs the approrpiate size for the arrays
	 * containing the results of the generalized eigen problem.
	 */
	void alloc(const int & nbrOFRows);

	/**
	 * This method frees the memory allocated.
	 */
	void free(void);

public:
	
	/**
	 * Destructor of the class.
	 */
	~generalizedEigenSolver(void);


	/**
	 * Method that clears the solver.
	 */
	void clear(void);
	
	/**
	 * Method to know if the solver is empty.
	 */
	bool empty(void) const;
	
	/**
	 * This method is use to acces the size of the problem that is solved.
	 *
	 * @return the number of eigen values (and eigen vectors) that are stored for the result.
	 */
	const int & getNbrOfEigenValues(void) const;

	/**
	 * This method is used to acces the eigen value.
	 */
	const double & getEigenValue(const int & item) const;

		/**
	 * This method is used to acces the eigen value.
	 */
	const containor<double> & getEigenValues(void) const;

	/**
	 * This method acces the eigen vectors.
	 */
	const vector & getEigenVector(const int & item) const;

	/**
	 * This method acces the eigen vectors.
	 */
	const containor<vector> & getEigenVectors(void) const;

	/**
	 * Méthode pour résoudre le problème.
	 */
	void solve(const matrixSymetric & A, const matrixSymetric & B);
};

#endif