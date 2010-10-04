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
#ifndef _GAUSSIAN_
#define _GAUSSIAN_

#include <containor.h>

class gaussianPolynome;

/**
 * Classe pour manipuler des fonctions gaussiennes.
 *
 * Cette classe permet de manipuler l'ensemble des fontctions de
 * la forme 
 *
 * \f C.exp(-\alpha(x-x)^2) \f
 *
 * Dans ce contexte C est le "coefficient" de la gaussienne,
 * x son "centre", et alpha son "exposant".
 *
 * Il est important de remarquer que le produit de deux gaussiennes 
 * est une gaussienne. Ceci permet de munir la classe d'un operateur
 * * qui effectue ce produit.
 */ 
class gaussian
{
private:
	
	/**
	 * Le centre de la gaussienne.
	 */
	double Center;
	
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
		* Méthode pour faire une copie de gaussiennes.
		*/
		void copy(const gaussian & g);
	
public:
		
		/**
		* Constructeur par défaut.
		 */
		gaussian(void);
	
	/**
	 * Constructeur avec une caractérisation de l'objet.
	 *
	 * @param coefficient le coefficient qui se trouve devant la gaussienne.
	 *
	 * @param exposant l'exposant de la gaussienne.
	 *
	 * @param center le centre de la gaussienne.
	 */
	gaussian(const double & coefficient ,const double & exponent ,const double & center);
	
	/**
		* Constructeur de copie.
	 *
	 * @param g la gaussienne avec laquelle on souahite initialiser l'objet crée.
	 */
	gaussian(const gaussian & g);
	
	/**
	 * Destructeur.
	 */
	virtual ~gaussian(void);
	
	/**
		* Méthode pour comparer deux gaussiennes :
	 * @TO DEFINE.
	 */
	static double compare(const gaussian & ga , const gaussian & gb);
	
	/**
		* Méthode pour calculer le centre de la gaussienne résultat du produit de deux gaussiennes.
	 *
	 * @return le centre de la gaussienne obtenus par la multiplication.
	 */
	static double computeCenter(const gaussian & ga , const gaussian & gb);
	
	/**
		* Méthode pour calculer le coefficient de la gaussienne issu du produit de deux gaussiennes.
	 * 
	 * @return le coefficient de la gaussienne obtenus par la multiplication.
	 */
	static double computeCoefficient(const gaussian & ga , const gaussian & gb);
	
	/**
		* Méthode pour calculer l'exposant de la gaussienne résultat du produit de deux gaussiennes.
	 *
	 * @return l'exposant de la gaussienne obtenus par multiplication
	 */
	static double computeExponent(const gaussian & ga , const gaussian & gb);
	
	/**
		* Méthode pour calculer la dérivée d'une gaussienne.
	 *
	 * @return le polynome gaussien qui est la dérivé de la gaussienne.
	 */
	gaussianPolynome derivate(void) const;
	
	/**
		* Méthode pour faire la division d'une gaussienne par un scalaire.
	 *
	 * Cette méthode permet d'effectuer la division d'une gaussienne
	 * et d'un scalaire.
	 * 
	 * @param g la gaussienne.
	 *
	 * @param scalar le scalaire par lequel on souahite diviser la gaussienne.
	 */
	void divide(const gaussian & g , const double & scalar);
	
	/**
	 * Méthode pour faire la division unaire d'une gaussienne par un scalaire.
	 *
	 * Cette méthode permet d'effectuer la division d'une gaussienne
	 * et d'un scalaire.
	 * 
	 * @param g la gaussienne.
	 *
	 * @param scalar le scalaire par lequel on souahite diviser la gaussienne.
	 */
	void divide(const double & scalar);
	
	/**
		* Méthode pour connaitre le centre de la gaussienne.
	 *
	 * @return le centre de la gaussienne.
	 */
	const double & getCenter(void) const;
	
	/**
		* Méthode pour connaitre le coefficient de la gaussienne.
	 *
	 * @return le coefficient de la gaussienne.
	 */
	const double & getCoefficient(void) const;
	
	/**
		* Méthode pour connaitre l'exposant de la gaussienne.
	 *
	 * @return l'exposant de la gaussienne.
	 */
	const double & getExponent(void) const;
	
	/**
		* Méthode qui calcule l'intégrale sur R de la gaussienne.
	 */
	double integral(void) const;
	
	/**
		* Méthode qui calcule le recouvrment pour la gaussienne et
	 * le produit de la gaussienne et de monomes de même centre.
	 */
	void integrals(const int & d_max , containor<double> & overlaps) const; 
	
	/**
		* Méthode qui calcule le recouvrment pour la gaussienne et
	 * le produit de la gaussienne et de monomes de même centre.
	 *
	 * Comme le calcul donne toujours 0 pour les degrés impairs
	 * cette méthode ne calcul que les degré pairs et renvoie un
	 * tableau de taille (d_max / 2  + 1) qui contient à l'emplacement
	 * la valeur pour un degré 2*i.
	 */
	void integrals_odd(const int & d_max , containor<double> & overlaps) const; 
	
	/**
		* Méthode pour faire la multiplication ed deux gaussienne.
	 *
	 * Cette méthode permet d'effectuer le produit de deux gaussienne
	 * et stoque le résultat dans l'objet appelant.
	 * 
	 * @param ga le première gaussienne.
	 *
	 * @param gb la seconde gaussienne.
	 */
	void multiply(const gaussian & ga , const gaussian & gb);
	
	/**
		* Méthode pour faire la unaire multiplication ed deux gaussienne.
	 *
	 * Cette méthode permet d'effectuer le produit de deux gaussiennes
	 * et stoque le résultat dans l'objet appelant.
	 *
	 * @param g la seconde gaussienne avec laquelle on souhaite multiplier l'objet appelant.
	 */
	void multiply(const gaussian & g);
	
	/**
		* Multiplication par un scalaire bianire.
	 *
	 * Méthode pour faire la multiplication de la gaussienne
   * avec un scalaire. Le résultat de l'opération est stoqué
   * dans l'objet appelant.
	 *
	 * @param g la gaussienne.
	 *
	 * @param scalar le scalaire.
	 */
	void multiply(const gaussian & g , const double & scalar);
	
	/**
		* Multiplication avec un scalaire unaire.
	 *
	 * Méthode pour faire la multiplication de l'objet avec un scalaire.
	 * 	 
	 * @param scalar le scalaire.
	 */
	void multiply(const double & scalar);
	
	/**
		* Opérateur d'affectation.
	 *
	 * @param g la gaussienne que l'on souahite recopier dans l'objet.
	 */
	gaussian & operator=(const gaussian & g);
	
	/**
		* Opérateur pour faire la multiplication de deux gaussiennes.
	 *
	 * @param g la gaussienne avec laquelle on souhaite multiplier l'objet.
	 *
	 * @return la valeur du produit.
	 */
	gaussian operator * (const gaussian & g) const;
	
	/**
		* Opérateur pour faire la multiplication d'une gaussienne et d'un scalaire.
	 *
	 * @param scalar le scalaire avec lequel on souhaite multiplier l'objet.
	 *
	 * @return la valeur du produit.
	 */
	gaussian operator * (const double & scalar) const;
	
	/**
		* Opérateur pour faire la multiplication unaire de deux gaussiennes.
	 *
	 * @param g la gaussienne avec laquelle on souhaite multiplier l'objet.
	 *
	 * @return la valeur du produit.
	 */
	gaussian & operator *= (const gaussian & g);
	
	/**
		* Opérateur pour faire la multiplication unaire d'une gaussienne et d'un scalaire.
	 *
	 * @param scalar le scalaire avec lequel on souhaite multiplier l'objet.
	 *
	 * @return la valeur du produit.	 
	 */
	gaussian & operator *= (const double & scalar);
	
	/**
	 * Opérateur pour faire la division d'une gaussienne et d'un scalaire.
	 *
	 * @param scalar le scalaire avec lequel on souhaite diviser l'objet.
	 *
	 * @return la gaussienne produit.
	 *
	 * @warning il vaut mieux que le scalaire soit non nul.
	 */
	gaussian operator / (const double & scalar) const;
	
	/**
	 * Opérateur pour faire la division unaire d'une gaussienne et d'un scalaire.
	 *
	 * @param scalar le scalaire avec lequel on souhaite diviser l'objet.
	 *
	 * @return une référence vers la gaussienne produit.
	 *
	 * @warning il vaut mieux que le scalaire soit non nul.
	 */
	gaussian & operator /= (const double & scalar);
	
	/**
		* Méthode qui calcule la valeur en x de la primitive de
	 * la gaussienne qui s'annule en 0.
	 */
	double primitive(const double & x) const;
	
	/**
		* Méthode qui calcule le recouvrment pour la gaussienne et
	 * le produit de la gaussienne et de monomes de même centre.
	 */
	void primitives(const double & x , const int & d_max , containor<double> & primitives) const; 
	
	/**
		* Méthode pour fixer le centre de la gaussienne.
	 */
	void setCenter(const double & center);
	
	/**
		* Méthode pour fixer le coefficient de la gaussienne.
	 */
	void setCoefficient(const double & coefficient);
	
	/**
		* Méthode pour fixer l'exposant de la gaussienne.
	 *
	 * @warning il faut que l'argument passé à cette fonction
	 * soit strictement positif.
	 */
	void setExponent(const double & exponent);		
	
	/**
		* Méthode pour écrire une gaussienne dans un flux.
	 */
	void write(ostream & outStream) const;
};

/**
 * Opérateur externe pour la multiplication d'un scalaire et d'une gaussienne.
 *
 * @param scalar le scalaire.
 *
 * @param g la gaussienne=
 *
 * @return une gaussienne qui contient le résultat du produit de la
 * gaussienne avec le scalaire.
 */
 extern gaussian operator* (const double & scalar , const gaussian & g);

/**
 * Opérateur externe d'écriture dans un flux.
 *
 * @param outStream le flux dans lequel on souhaite écrire la gaussienne.
 *
 * @param g la gaussienne à écrire.
 */
extern ostream & operator<<(ostream & outStream , const gaussian & g);


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
#include "gaussianPolynome.h"

#endif
