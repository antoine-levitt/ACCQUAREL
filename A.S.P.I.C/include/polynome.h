/* 
* The polynome library of the A.S.P.I.C. 
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
#ifndef _POLYNOME_
#define _POLYNOME_
#include "polynomeBase.h"


/**
 * Classe pour la manipulation des polynomes standarts.
 * 
 * Usuellement , un polynome peut s'ecrire comme p(x) = sum_{i}{a_i * x^i}. Le propos de cette
 * classe est de permettre la manipulation de ces quantités. 
 * Cette classe définie pour ces objets les principaux opérateurs qu'elle peut utiliser.
 *
 * une addition interne :  si p et q sont deux polynomes , il est possible d'écrire p+q
 * dans un programme pour effectuer la sommation de ces deux polynomes.
 * @see polynome::operator+ (const polynome &) const;
 * @see polynome::operator+= (const polynome &);
 * 
 * une soustraction interne : si p et q sont deux polynomes , il est possible d'écrire p-q
 * dans un programme pour effectuer la soustraction de ces deux polynomes.
 * @see polynome::operator- (const polynome &) const;
 * @see polynome::operator-= (const polynome &);
 *
 * une multiplication interne.
 * @see polynome::operator* (const polynome &) const;
 * @see polynome::operator*= (const polynome &);
 *
 * une multiplication avec les scalaire.
 * @see polynome::operator* (const double &) const;
 * @see polynome::operator-= (const double &);
 *
 * une division avec des scalaires.
 * @see polynome::operator/ (const double &) const;
 * @see polynome::operator/= (const double &);
 *
 * La classe propose aussi une méthode pour effectuer l'opération de dérivation 
 * d'un polynome par rapport à sa variable.
 * @see polynome::derviate().
 *
 *
 */
class polynome : public polynomeBase<double>
{
private:

protected:

	/**
	 * Methode pour afficher les fonction de base.
	 *
	 * Comme nous travaillons ici avec les polynômes usuel la 
	 * fonction de base pour un degré d, s'écriit "x^d".
	 *
	 * @param i le degré de la fonction de base qui l'on cherche à convertir en
	 * chaine de charactères.
	 *
	 * @return la chaine de charactère qui repésente la fonction de base, ie "1" lorsque
	 * le degré est nul et "x^d" lorsque le degré demandé est supérieur ou égal à 1.
	 */
	virtual const string  getBaseString(const int & i) const;

public:

	/**
	 * Constructeur par defaut.
	 */
	polynome(void);

	/**
	 * Constructeur de copie.
	 *
	 * @param poly le polynome avec lequel on souhaite initialiser l'objet.
	 */
	polynome(const polynome & poly);

	/**
	 * Construcuteur de polynome constant.
	 *
	 * Ce constructeur permet de créer un polynome constant dont la
	 * valeur est passée en argument.
	 */
	polynome(const double & value);

	/**
	* Destructeur.
	*/
	virtual ~polynome(void);

	/**
	 * Methode qui multiplie un polynome par (x-root).
	 * Autrement dit cela revient à ajouter root comme zéro du polynome.
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @return le polynome produit de (x-roo) et de l'objet.
	 */
	polynome centerMonomeMultiply(const double & root)const;

	/**
	 * Methode qui multiplie un polynome par (x-root)^degree.
	 * Autrement dit cela revient à ajouter root comme zéro du polynome
	 * avec un ordre degree.
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @param degree l'ordre avec lequel on souhaite ajouter la racine au polynome. 
	 *
	 * @return le polynome produit de (x-roo) et de l'objet.
	 */
	 polynome centerMonomeMultiply(const double & root , const int & degree) const;

	/**	
	 * Methode pour dériver le polynome.
	 *
	 * @return le polynome issu de la dérivation de l'objet.
	 */
	 polynome derivate(void) const;

	/**
	 * Methode qui multiplie un polynome par x^degree.
	 *
	 * @param degree le degree du monome avec lequel souhaite multiplier le polynome.
	 *
	 * @return le polynome qui contient le produit de x^d et de l'objet.
	 */
	polynome monomeMultiply(const int & degree=1) const;

	/**
	 * Operateur pour l'assignement.
	 * 
	 * @param p le polynome que l'on souhaite copier dans
	 * l'objet appelant.
	 *
	 * @return une référence vers le nouveau polynome.
	 */
	polynome & operator= (const polynome & p);

	/**
	 * Operateur pour faire l'addition de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite ajouter à l'objet.
	 *
	 * @return un polynome qui contient la somme demandée.
	 */
	polynome operator+ (const polynome & p) const;

	/**
	 * Opérateur pour faire une addition unaire de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite ajouter à l'objet.
	 *
	 * @return une référence vers le polynome qui contient la somme demandée.
	 */
	polynome & operator+= (const polynome & p);	
	
	/**
	 * Operateur pour faire la soustraction de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite soustraire ajouter à l'objet.
	 *
	 * @return un polynome qui contient la différence demandée.
	 */
	polynome operator- (const polynome & p) const;

	/**
	 * Operateur pour faire la soustraction unaire de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite soustraire à l'objet.
	 *
	 * @return une référence vers le polynome qui contient la différence demandée.
	 */
	polynome & operator-= (const polynome & p);

	/**
	 * Operateur pour faire la multiplication de deux polynomes.
	 * 
	 * @param p le polynome avec lequel on souhaite multiplier l'objet.
	 *
	 * @return un polynome qui contient le produit demandé.
	 */
	polynome operator* (const polynome & p) const;

	/**
	 * Operateur pour faire la multiplication d'un polynome et d'un scalaire.
	 * 
	 * @param scalar le scalaire avec lequel on souhaite multiplier l'objet.
	 *
	 * @return un polynome qui contient le produit demandé.
	 */
	polynome operator* (const double & scalar) const;

	/**
	 * Operateur pour faire une multiplication unaire de deux polynomes.
	 * 
	 * @param p le polynome avec lequel on souhaite multiplier l'objet.
	 *
	 * @return une référence vers le polynome qui contient le produit demandé.
	 */
	polynome & operator*= (const polynome & p);

	/**
	 * Operateur pour faire une multiplication unaire d'un polynome et d'un scalaire.
	 * 
	 * @param scalar le scalaire avec lequel on souhaite multiplier l'objet.
	 *
	 * @return une référence vers le polynome qui contient le produit demandé.
	 */
	polynome & operator *= (const double & scalar);

	/**
	 * Operateur pour faire la disision d'un polynome et d'un scalaire.	
	 *
	 * @param scalar le scalaire par lequel on souhaite diviser l'objet.
	 *
	 * @return un polynome qui contient le quotient demandé.
	 */
	polynome operator/ (const double & scalar) const;

	/**
	 * Operateur pour faire une division unaire d'un polynome et d'un sclaire.
	 *
	 * @param scalar le scalaire par lequel on souhaite diviser l'objet.
	 *
	 * @return une référence vers le polynome qui contient le quotient demandé.
	 */
	polynome & operator/= (const double & scalar);
	
	/**
	 * Méthode pour changer la varialble linéaire.
	 *
	 * @param scalar le facteur du changement de variable.
	 *
	 * @return la forme du polynome après le changement de variable.
	 */
	polynome scale(const double & scalar);

	/**
	 * Méthode pour donner au polynome une valeur constante.
	 *
	 * @param value la valeur que l'on souhaite donner au polynome.
	 */
	void setPolynomeCoefficients(const double & value);

	/**
	 * Méthode pour recopier la valeur d'un polynome.
	 *
	 * @param poly la valeur que l'on souhaite recopier.
	 *
	 * @see operator= (const polynome & poly).
	 */
	void setPolynomeCoefficients(const polynome & poly);

};

/**
 * Opérateur de multiplication d'un scalaire et d'un polynome.
 *
 * @param scalar un scalaire.
 *
 * @param poly un polynome
 *
 * @return le polynome produit du scalaire et du polynome.
 */
extern polynome operator*(const double & scalar , const polynome & poly);

#endif
