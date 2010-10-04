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
#ifndef _POLYNOME_COMPLEX_
#define _POLYNOME__COMPLEX_

#include "polynomeBase.h"
#include <complex.h>


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
 * @see polynome::operator* (const complex &) const;
 * @see polynome::operator-= (const complex &);
 *
 * une division avec des scalaires.
 * @see polynome::operator/ (const complex &) const;
 * @see polynome::operator/= (const complex &);
 *
 * La classe propose aussi une méthode pour effectuer l'opération de dérivation 
 * d'un polynome par rapport à sa variable.
 * @see polynome::derviate().
 *
 *
 */
class complexPolynome : public polynomeBase<complex>
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
	complexPolynome(void);

	/**
	 * Constructeur de copie.
	 *
	 * @param poly le polynome avec lequel on souhaite initialiser l'objet.
	 */
	complexPolynome(const complexPolynome & poly);

	/**
	 * Construcuteur de polynome constant.
	 *
	 * Ce constructeur permet de créer un polynome constant dont la
	 * valeur est passée en argument.
	 */
	complexPolynome(const complex & value);

	/**
	* Destructeur.
	*/
	virtual ~complexPolynome(void);

	/**
	 * Methode qui multiplie un polynome par (x-root).
	 * Autrement dit cela revient à ajouter root comme zéro du polynome.
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @return le polynome produit de (x-roo) et de l'objet.
	 */
	complexPolynome centerMonomeMultiply(const complex & root)const;

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
	 complexPolynome centerMonomeMultiply(const complex & root , const int & degree) const;

	/**	
	 * Methode pour dériver le polynome.
	 *
	 * @return le polynome issu de la dérivation de l'objet.
	 */
	 complexPolynome derivate(void) const;

	/**
	 * Methode qui multiplie un polynome par x^degree.
	 *
	 * @param degree le degree du monome avec lequel souhaite multiplier le polynome.
	 *
	 * @return le polynome qui contient le produit de x^d et de l'objet.
	 */
	complexPolynome monomeMultiply(const int & degree=1) const;

	/**
	 * Operateur pour l'assignement.
	 * 
	 * @param p le polynome que l'on souhaite copier dans
	 * l'objet appelant.
	 *
	 * @return une référence vers le nouveau polynome.
	 */
	complexPolynome & operator= (const complexPolynome & p);

	/**
	 * Operateur pour faire l'addition de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite ajouter à l'objet.
	 *
	 * @return un polynome qui contient la somme demandée.
	 */
	complexPolynome operator+ (const complexPolynome & p) const;

	/**
	 * Opérateur pour faire une addition unaire de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite ajouter à l'objet.
	 *
	 * @return une référence vers le polynome qui contient la somme demandée.
	 */
	complexPolynome & operator+= (const complexPolynome & p);	
	
	/**
	 * Operateur pour faire la soustraction de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite soustraire ajouter à l'objet.
	 *
	 * @return un polynome qui contient la différence demandée.
	 */
	complexPolynome operator- (const complexPolynome & p) const;

	/**
	 * Operateur pour faire la soustraction unaire de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite soustraire à l'objet.
	 *
	 * @return une référence vers le polynome qui contient la différence demandée.
	 */
	complexPolynome & operator-= (const complexPolynome & p);

	/**
	 * Operateur pour faire la multiplication de deux polynomes.
	 * 
	 * @param p le polynome avec lequel on souhaite multiplier l'objet.
	 *
	 * @return un polynome qui contient le produit demandé.
	 */
	complexPolynome operator* (const complexPolynome & p) const;

	/**
	 * Operateur pour faire la multiplication d'un polynome et d'un scalaire.
	 * 
	 * @param scalar le scalaire avec lequel on souhaite multiplier l'objet.
	 *
	 * @return un polynome qui contient le produit demandé.
	 */
	complexPolynome operator* (const complex & scalar) const;

	/**
	 * Operateur pour faire une multiplication unaire de deux complexPolynomes.
	 * 
	 * @param p le polynome avec lequel on souhaite multiplier l'objet.
	 *
	 * @return une référence vers le polynome qui contient le produit demandé.
	 */
	complexPolynome & operator*= (const complexPolynome & p);

	/**
	 * Operateur pour faire une multiplication unaire d'un polynome et d'un scalaire.
	 * 
	 * @param scalar le scalaire avec lequel on souhaite multiplier l'objet.
	 *
	 * @return une référence vers le polynome qui contient le produit demandé.
	 */
	complexPolynome & operator *= (const complex & scalar);

	/**
	 * Operateur pour faire la disision d'un polynome et d'un scalaire.	
	 *
	 * @param scalar le scalaire par lequel on souhaite diviser l'objet.
	 *
	 * @return un polynome qui contient le quotient demandé.
	 */
	complexPolynome operator/ (const complex & scalar) const;

	/**
	 * Operateur pour faire une division unaire d'un polynome et d'un sclaire.
	 *
	 * @param scalar le scalaire par lequel on souhaite diviser l'objet.
	 *
	 * @return une référence vers le polynome qui contient le quotient demandé.
	 */
	complexPolynome & operator/= (const complex & scalar);
	
	/**
	 * Méthode pour changer la varialble linéaire.
	 *
	 * @param scalar le facteur du changement de variable.
	 *
	 * @return la forme du polynome après le changement de variable.
	 */
	complexPolynome scale(const complex & scalar);
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
extern complexPolynome operator*(const complex & scalar , const complexPolynome & poly);

#endif
