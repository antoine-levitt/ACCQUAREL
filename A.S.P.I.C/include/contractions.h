/* 
* The chemics library of A.S.P.I.C. 
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
#ifndef _CONTRACTION_
#define _CONTRACTION_

#include <containor.h>

/**
 * Cette classe permet de gÃ©rer les contractions.
 */
class contractions 
{
private:
	
	/*
   * array that contains the coefficients of the contraction.
	 */
	containor<double> Coefficients;

	/**
	 * array that contains the exponents of the contraction.
	 */
	containor<double> Exponents;

//protected:
/**
 * On rend cette partie de la classe "contractions" publique afin de pouvoir interfacer le code A.S.P.I.C.
 * avec un code écrit en fortran95/2003.
 * (modification de G. Legendre (CEREMADE - université Paris-Dauphine) le 12/07/2007)
 **/

public:
	
	/**
	 * Method COPY.
	 */
	void copy(const contractions & c);

	/**
	 * Method to copy directly the array of coefficients.
	 *
	 * @warning this method does not performs memory allocation.
	 * the number of contractions must be set to the corrct value
	 * before the call of this method.
	 */
	void setCoefficients(const containor<double> & coefficients);

	/**
	 * Method to copy directly the array of exponents.
	 *
	 * @warning this method does not performs memory allocation.
	 * the number of contractions must be set to the corrct value
	 * before the call of this method.
	 *
	 * @warning the array of exponents pass as argument must contains
	 * only strictly positive values.
	 */
	void setExponents(const containor<double> & coefficients);

//public:

	/**
	 * default constructor.
	 */
	contractions(void);

	/**
	 * Copy constructor.
	 */
	contractions(const contractions & c);
	
	/**
	 * Destructor.
	 */
	virtual ~contractions(void);

	/**
	 * Method that clears the object.
	 */
	void clear(void);
	

	/**
	 * Méthode pour savoir si l'objet contractions est vide.
	 */
	bool empty(void) const;

	/**
   * Method GET for a contraction coefficient.
	 *
	 * @param item the number of the contraction for which
	 * the exponent will be returned.
	 *
	 * @warning the item parameter must be non negative
	 * and strictly inferior than the value returned by the method
	 * getNbrOfContractions().
	 */
	const double & getCoefficient(const int & item) const;
	
	/**
	 * Methode GET for all coefficents.
	 */
	const containor<double> & getCoefficients(void) const;

	/**
	 * Method GET for the exponent of the contraction..
	 *
	 * @param item the number of the contraction for which
	 * the exponent will be returned.
	 *
	 * @warning the item parameter must be non negative
	 * and strictly inferior than the value returned by the method
	 * getNbrOfContractions().
	 */
	const double & getExponent(const int & item) const;
	
	/**
	 * Method GET for all exponents.
	 */
	const containor<double> & getExponents(void) const;

	/**
   * Methode GET for the number of contractions.
	 */
	int getNbrOfContractions(void) const;
	
	/**
	 * Operator =
	 */
	contractions & operator=(const contractions & c);
	
	/**
	 * Method SET for the coefficient.
	 * 
	 * @param item the number of the contraction for which
	 * the exponent will be returned.
	 *
	 * @warning the item parameter must be non negative
	 * and strictly inferior than the value returned by the method
	 * getNbrOfContractions().
	 */
	void setCoefficient(const int & item , const double & coefficient);

	/**
	 * Method SET for the all the contractions.
	 */
	void setContractions(const contractions & c);

	/**
	 * Method SET for an exponent of the contraction.
	 *
	 * @param item the number of the contraction for which
	 * the exponent will be returned.
	 *
	 * @warning the item parameter must be non negative
	 * and strictly inferior than the value returned by the method
	 * getNbrOfContractions().
	 */
	void setExponent(const int & item , const double & exponent);
	
	/**
   * Method SET for the number of contractions.
	 */
	void setNbrOfContractions(const int & nbrOfContractions);

	/**
	 * Method write.
	 */
	void write(ostream & out) const;
};

/**
 * Operator <<.
 */
ostream & operator<<(ostream & out, const contractions & c);

#endif
