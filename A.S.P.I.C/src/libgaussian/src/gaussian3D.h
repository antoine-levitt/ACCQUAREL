/* 
* The gaussian library of the A.S.P.I.C. 
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
#ifndef _GAUSSIAN_3D_
#define _GAUSSIAN_3D_

#include <containor3D.h>
#include <dpoint.h>

/**
 * On déclare ici la classe gaussianPolynome_3D.
 */
class gaussianPolynome3D;

class gaussian3D
{
private:
	/** 
	 * Le centre de la gaussienne.
	 */
	dpoint<3> Center;

	/**
	 * Le coefficient devant la gaussienne.
	 */
	double Coefficient;

	/**
	 * L'exposant de la gaussienne.
	 */
	double Exponent;
	
protected:

	/**
	 * Méthode pour copier les gaussiennes.
	 */
	void copy(const gaussian3D & g);	

public:

	/**
	 * Constructeur de la classe gaussienne.
	 */
	gaussian3D(void);
	
	/**
	 * Constrcuteur avec spécification des paramêtres.
	 */
	gaussian3D(const dpoint<3> & center , const double & coefficient , const double & exponent);

	/**
	 * Constructeur de copie.
	 */
	gaussian3D(const gaussian3D & g);

	/**
	 * Destructeur de la classe gaussienne.
	 */
	~gaussian3D(void);

	/**
	 * Méthode pour comparer deux gaussiennes :
	 * @TO DEFINE.
	 */
	static double compare(const gaussian3D & ga , const gaussian3D & gb);
	
	/**
	 * Méthode pour calculer le centre de la gaussienne
	 * produit de deux gaussiennes.
	 */
	static dpoint<3> computeCenter(const gaussian3D  & ga , const gaussian3D & gb);


	/**
	 * Méthode pour calculer le coefficient de la gaussienne
	 * produit de deux gaussiennes.
	 */
	static double computeCoefficient(const gaussian3D & ga , const gaussian3D & gb);

	/**
	 * Méthode pour calculer le coefficient de la gaussienne
	 * produit de deux gaussiennes.
	 */
	static double computeExponent(const gaussian3D & ga , const gaussian3D & gb);

	/**
	 * Méthode pour calculer la dérivée d'une gaussienne.
	 *
	 * Attention !! le parametre dim représente la direction
	 * de dérivation : 0 repésente la dérivation par rapport
	 * à x, 1 par rapport à y, 2 par rapport à z.
	 */
	gaussianPolynome3D derivate(int dim) const;

	/**
	 * Method to eval the value of the gaussian 
	 * in the point x.
	 */
	double eval(const dpoint<3> & position) const;
	
	/**
	 * Méthode pour connaitre le centre.
	 */
	const dpoint<3> & getCenter(void) const;
	
	/**
	 * Méthode pour connaitre le coefficient.
	 */
	const double & getCoefficient(void) const;

	/**
	 * Méthode pour connaitre l'exposant.
	 */ 
	const double & getExponent(void) const;

	/**
	 * Méthode pour calculer l'intégrale sur R de la gaussienne.
	 *
	 * \f[ 
	 * \int_{R}{C \exp(-\alpha (X-X_0)^2 )}
	 * \f]
	 */
	double integral(void) const;
	
	/**
	 * Méthode pour calculer l'intégrale sur R de la gaussienne.
	 */
	void integrals(const ipoint<3> & deg_max , containor3D<double> & overlaps) const;
	
	/**
	 * Méthode pour faire la multiplication.
	 */
	void multiply(const gaussian3D & ga , const gaussian3D & gb);

	/**
	 * Méthode pour faire la multiplication.
	 */
	void multiply(const gaussian3D & g);
	
	/**
	 * Méthode pour faire la multiplication.
	 */
	void multiply(const gaussian3D & g , const double & scalar);

	/**
	 * Méthode pour faire la multiplication.
	 */
	void multiply(const double & scalar);
	
	/**
	 * Multiplication sous forme d'opérateur.
	 */
	gaussian3D operator*(const gaussian3D & g) const;
	
	/**
	 * Multiplication sous forme d'opérateur.
	 */
	gaussian3D & operator*=(const gaussian3D & g);
	
	/**
	 * Multiplication sous forme d'opérateur.
	 */
	gaussian3D operator*(const double & scalar) const;
	
	/**
	 * Multiplication sous forme d'opérateur.
	 */
	gaussian3D & operator*=(const double & scalar);
	
	/**
	 * Méthode pour fixer le centre.
	 */
	void setCenter(const dpoint<3> & center);

	/**
	 * Méthode pour fixer le centre.
	 */
	void setCenter(const double & center_x, const double & center_y , const double & center_z);

	/**
	 * Méthode pour fixer le coefficient.
	 */
	void setCoefficient(const double & coefficient);

	/**
	 * Méthode pour fixer l'exposant.
	 */
	void setExponent(const double & exponent);

	/**
	 * Méthode pour écrire.
	 */
	void write(ostream & outStream) const;
	
};

/**
 * Opérateur externe de multiplication par un scalaire.
 */
extern gaussian3D operator*(const double & scalar , const gaussian3D & g);


/**
 * Opérateur externe d'écriture dans un flux.
 */
extern ostream & operator<<(ostream & outStream , const gaussian3D & g);


/////////////////////////////////////////////////////////////////////////////////////
//  !! DO NOT MOVE !!
//
// On inclut les déclaration de classe gaussianPolynome3D pour que l'utilisateur 
// puisse utiliser la méthode dérivate aisément (pas vraiment d'include à rajouter
// dans le bon ordre dans le programe appelant).
//
// Mais il faut mettre cet include après les déclaration de la classe gaussian3D 
// parce que la compilation de la classe gaussianPolynome3D nécessite 
// la connaissance de la classe gaussienne
////////////////////////////////////////////////////////////////////////////////////
#include "gaussianPolynome3D.h"


#endif

