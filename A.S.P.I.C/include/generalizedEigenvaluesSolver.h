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

#include "vector.h"
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
class generalizedEigenvaluesSolver
{
private:

protected:

public:
	
	/**
	 * Méthode pour accéder aux valeurs propres.
	 */
	const double & getEigenvalue(const int & item) const;
	
	/**
	 * Méthode pour accéder aux vecteurs propres.
	 */
	const vector & getEigenvector(const int & item) const;
	
	/**
	 * Méthode pour résoudre le problème.
	 */
	void solve(const matrixSymetric & A, const matrixSymetric & B);
};