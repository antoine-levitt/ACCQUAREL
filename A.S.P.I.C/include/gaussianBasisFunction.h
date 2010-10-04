/* 
 * The chemics library of A.S.P.I.C. 
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
#ifndef	_GAUSSIAN_BASIS_FUNCTION_
#define	_GAUSSIAN_BASIS_FUNCTION_
#include <containor.h>
#include "contractions.h"
#include <dpoint.h>
#include <gaussianPolynome3D.h>
#include <ipoint.h>

/**
*	Classe qui repésente les fonctions de	base 
*	gaussienne.	
* 
* gaussianBasisFunction is
* the product of a momome in the three variables $x-x_c,y-y_c,z-z_c$ by a linear combination
* of gaussian functions centered at $(x_c,y_c,z_c)$
* thus a gaussianBasisFunction  is completely defined by
* its center $(x_c,y_c,z_c)$
* a containor containing the exponents of the gaussian functions in the linear combination
* a containor containing the coefficients  of the gaussian functions in the linear combination
* the degrees of the polynomial function, i.e. three integers
*/

class	gaussianBasisFunction	:	public contractions
{
private:

	/**	 
	 *	The center of the basis function.
	 */
	dpoint<3> Center;

	/**	
	 * the degree of the monome in the basis function
	 */	
	ipoint<3>	MonomeDegree;	

protected:	
	
		/**
	   * Method that performs the copy of a basis function.
	   */	
	void copy(const	gaussianBasisFunction	&	gbf);

public:	
	/**
	 * Default constructor.
	 */	
	gaussianBasisFunction(void);
	
	/**
	 * Copy constructor.
	 */	
	gaussianBasisFunction(const	gaussianBasisFunction	&	gbf);
	
	/**
	 * Destructor.
	 */	
	virtual	~gaussianBasisFunction(void);
	
	/**
	 * Method clear.
	 */
	void clear(void);
	
	/**
	 * Method to eval the value of the basis function 
	 * in the point x.
	 */
	double eval(const dpoint<3> & position) const;
	
	/**
	 * Method GET for the center.
	 */	
	const	dpoint<3>	&	getCenter(void)	const;
	
	/**
	 * Method GET to access the component of the 
	 * basis function as a gaussian polynome.
	 * 
	 * @warning the parameter item must be in the range of
	 * the contractions sizes.
	 */	
	gaussianPolynome3D getGaussianPolynome3D(const int & item) const;
	
	/**
	 * Method GET for the degree of the polynome.
	 */	
	const	ipoint<3>	&	getMonomeDegree(void) const;	
	
	/**
	 * Method SET for the center of the basis function.
	 */	
	void setCenter(const dpoint<3> & center);
	
	/**
	 * Methos SET for the degree of the polynome.
	 *
	 * @warning all the components of the monomeDegree object
	 * must be non negative integers.
	 */	
	void setMonomeDegree(const	ipoint<3>	&	monomeDegree);	
	
	/**	 
	 * Method SET for the degree of the polynome.
	 *
	 * @warning all the components of the monomeDegree object
	 * must be non negative integers.
	 */	
	void setMonomeDegree(int degree_x	,	int	degree_y	, int degree_z);	

	/**
	 * Method write.
	 */
	void write(ostream & out) const;

};

/**
 * Extern operator <<.
 */
inline ostream & operator<< (ostream & out , const gaussianBasisFunction & gbf)
{
	gbf.write(out);
	return out; 
}

#endif

