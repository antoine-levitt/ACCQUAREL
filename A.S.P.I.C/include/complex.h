/* 
 * The complex library of the A.S.P.I.C. 
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
 *
 * CERMICS, ENPC, hereby disclaims all copyright interest in
 * the library `complex' (a library of complex for A.S.P.I.C.) written
 * by François Lodier.
 */
#ifndef _COMPLEX_
#define _COMPLEX_

#include <iostream>
#include <math.h>

using namespace std;

/**
 * Classe pour la manipulation des nombres complexes.
 */
class complex 
{private:	
	/**
	 * La partie réelle du nombre complexe.
	 */
	double Real;
	
	/**
	 * La partie imaginaire du nombre complexe.
   */
	double Im;protected:
	
	/**
	 * Méthode pour copier un nombre complexe. 
	 */
	void copy(const complex & c);public:
	
	/**
	 * Constructeur par défaut.	
   */
	complex(void);	
	
	/**
	 * Constructeur à partir d'un nombre réel.
	 *
	 * @param real le nombre réel qui est transformé en complexe
	 * par l'appel de ce contructeur.
	 */
	complex(double real);	
	
	/**
	 * constructeur avec un reel et un imaginaire.
	 *
	 * @param real la partie réelle du nombre complexe.
	 *
	 * @param im la partie imaginaire du nombre complexe.
	 */
	complex(double real , double im);	
	
	/**
	 * constrcteur de copie.
	 *
	 * @param c le nombre complexe à copier dans l'objet construit.
	 */
	complex(const complex & c);	
	
	/**
	 * Destructeur.
	 */
	~complex(void);	
	
	/**
	 * Méthode constante pour accéder à la partie réelle.
	 *
	 * @return une référence constante vers partie réelle de l'objet.
	 */
	const double & real(void) const;	
	
	/**
	 * Méthode pour accéder à la partie réelle.
	 *
	 * @return une référence vers la partie réelle de l'objet.
	 */
	double & real(void);	
	
	/**
	 * Méthode constante pour accéder à la partie imaginaire.
	 *
	 * @return une référence constante vers partie imaginaire de l'objet.
	 */
	const double & im(void) const;	
	
	/**
	 * Méthode pour accéder à la partie imaginaire.
	 *
	 * @return une référence vers partie imaginaire de l'objet.
	 */
	double & im(void);
	
	/**
	 * Méthode pour savoir si un nombre complexe est un imaginaire pur.
	 *
	 * @return true lorsque la partie reele est nulle et false sinon.
	 */
	bool is_pure_im(void) const;	
	
	/**
	 * Méthode pour savoir si l'objet est un nombre est un reel.
	 *
	 * @return true lorsque la partie réelle est nulle, false sinon.
	 */
	bool is_real(void) const;	
	
	/**
	 * Méthode pour additionner deux complexes.
	 *
	 * @param c le nombre complex que l'on souhaite ajouter à l'objet appelant.
	 *
	 * @return un complex qui contient le résultat de l'addition de l'objet courrant 
	 * et de c.
	 */
	complex operator+(const complex & c) const;	
	
	/**
	 * Méthode pour additionner un complexe et un double :
	 *
	 * @param real le nombre réel que l'on souhaite ajouter à l'objet appelant.
	 *
	 * @return un complex qui contient le résultat de l'addition de l'objet courrant 
	 * et de real.
	 */
	complex operator+(const double & real) const;	
	
	/**
	 * Méthode pour faire l'addition unaire de deux complexes.
	 * 
	 * @param c le nombre complex que l'on souhaite ajouter à l'objet appelant.
	 *
	 * @return une référence vers l'objet qui contient le résultat de l'addition 
	 * de c et de l'objet appelant.
	 */
	complex & operator+=(const complex & c);	
	
	/**
	 * Méthode pour faire l'addition unaire d'un complexe et d'un double.
	 * 
	 * @param real le nombre réel que l'on souhaite ajouter à l'objet appelant.
	 *
	 * @return une référence vers le complex qui contient le résultat de l'addition de l'objet courrant 
	 * et de real.
	 */
	complex & operator+=(const double & real);	
	
	/**
	 * Méthode pour soustraire deux complexes.
	 *
	 * @param c le nombre complex que l'on souhaite soustraire à l'objet appelant.
	 *
	 * @return un complex qui contient le résultat de la soustraction de l'objet courrant 
	 * et de c.
	 */
	complex operator-(const complex & c) const;	
	
	/**
	 * Méthode pour soustraire un complexe et un double :
	 *
	 * @param real le nombre réel que l'on souhaite soustraire à l'objet appelant.
	 *
	 * @return un complex qui contient le résultat de la soustraction de l'objet courrant 
	 * et de real.
	 */
	complex operator-(const double & real) const;	
	
	/**
	 * Méthode pour faire la soustraction unaire de deux complexes.
	 * 
	 * @param c le nombre complex que l'on souhaite soustraire à l'objet appelant.
	 *
	 * @return une référence vers l'objet qui contient le résultat de la soustraction 
	 * de c et de l'objet appelant.
	 */
	complex & operator-=(const complex & c);	
	
	/**
	 * Méthode pour faire la soustraction unaire d'un complexe et d'un double.
	 * 
	 * @param real le nombre réel que l'on souhaite soustraire à l'objet appelant.
	 *
	 * @return une référence vers le complex qui contient le résultat de la soustraction de l'objet courrant 
	 * et de real.
	 */
	complex & operator-=(const double & real); 
	
	/**
	 * Méthode pour multiplier deux complexes.
	 *
	 * @param c le nombre complex que l'on souhaite multiplier avec l'objet appelant.
	 *
	 * @return un complex qui contient le résultat de la multiplication de l'objet courrant 
	 * et de c.
	 */
	complex operator*(const complex & c) const;	
	
	/**
	 * Méthode pour multiplier un complexe et un double :
	 *
	 * @param real le nombre réel que l'on souhaite multiplier avec l'objet appelant.
	 *
	 * @return un complex qui contient le résultat de la multiplication de l'objet courrant 
	 * et de real.
	 */
	complex operator*(const double  & real) const;	
	
	/**
	 * Méthode pour faire la multiplication unaire de deux complexes.
	 * 
	 * @param c le nombre complex que l'on souhaite multiplier avec l'objet appelant.
	 *
	 * @return une référence vers l'objet qui contient le résultat de la multiplication 
	 * de c et de l'objet appelant.
	 */
	complex & operator*=(const complex & c);	
	
	/**
	 * Méthode pour faire la multiplication unaire d'un complexe et d'un double.
	 * 
	 * @param real le nombre réel que l'on souhaite multiplier avec l'objet appelant.
	 *
	 * @return une référence vers le complex qui contient le résultat de la multiplication de l'objet courrant 
	 * et de real.
	 */
	complex & operator*=(const double & real); 

	/**
	 * Méthode pour diviser deux complexes.
	 *
	 * @param c le nombre complex par lequel on souhaite diviser l'objet appelant.
	 *
	 * @return un complex qui contient le résultat de la division de l'objet courrant 
	 * par c.
	 */
	complex operator /(const complex & c) const;	
	
	/**
	 * Méthode pour diviser un complexe et un double :
	 *
	 * @param real le nombre réel par lequel souhaite diviser l'objet appelant.
	 *
	 * @return un complex qui contient le résultat de la division de l'objet courrant 
	 * par real.
	 */
	complex operator /(const double & real) const;	
	
	/**
	 * Méthode pour faire la division unaire de deux complexes.
	 * 
	 * @param c le nombre complex par lequel on souhaite diviser avec l'objet appelant.
	 *
	 * @return une référence vers l'objet qui contient le résultat de la division 
	 * de l'objet appelant par c.
	 */
	complex & operator/=(const complex & c);	
	
	/**
	 * Méthode pour faire la division unaire d'un complexe et d'un double.
	 * 
	 * @param real le nombre réel par lequel on souhaite diviser l'objet appelant.
	 *
	 * @return une référence vers le complex qui contient le résultat de la division de l'objet courrant 
	 * par real.
	 */
	complex & operator/=(const double & real);
	
	/**
	 * Opérateur d'affectation
	 * 
	 * @param c la valeur que l'on souhaite affecter à l'objet.
	 *
	 * @return une référence vers l'objet affecté de sa nouvelle valeur 
	 */
	complex & operator=(const complex & c);	
	
	/**
	 * Opérateur d'affectation
	 * 
	 * @param real la valeur que l'on souhaite affecter à l'objet.
	 *
	 * @return une référence vers l'objet affecté de sa nouvelle valeur 
	 */
	complex & operator=(const double & real);	
	
	/**
	 * Méthode pour calculer le module au carré d'un nombre complexe.
	 *
	 * @return la valeur du module au carré du nombre complexe.
	 */
	double sq_module(void) const;	
	
	/**
	 * Méthode pour calculer le module d'un nombre complexe.
	 *
	 * @return la valeur du module du nombre complexe.
	 */
	double module(void) const;

	/**
	 * Méthode pour calculer le conjugué d'un nombre complexe.
	 *
	 * @return le complex conjugué de l'objet.
	 */
	complex conjugate(void) const;	
	
	/**
	 * Opérateur de comparaison pour deux nombres complexes.
	 *
	 * @param c le nombre complexe avec lequel on souhaite comparer l'objet appelant.
	 *
	 * @return true lorsque c et l'objet sont identiques, flase sinon.
	 */
	bool operator==(const complex & c) const;	
	
	/**
	 * Opérateur de comparaison pour deux nombres complexes.
	 *
	 * @param c le nombre complexe avec lequel on souhaite comparer l'objet appelant.
	 *
	 * @return false lorsque c et l'objet sont identiques, true sinon.
	 */
	bool operator!=(const complex & c) const;	
	
	/**
	 * Opérateur de comparaison pour un nombre complexe et un reel.
	 *
	 * @param real le nombre réel avec lequel on souhaite comparer l'objet appelant.
	 *
	 * @return true lorsque r et l'objet sont identiques, false sinon.
	 */
	bool operator==(const double & real) const;	
	
	/**
	 * Opérateur de comparaison pour un nombre complexe et un reel.
	 *
	 * @param real le nombre réel avec lequel on souhaite comparer l'objet appelant.
	 *
	 * @return false lorsque r et l'objet sont identiques, true sinon.
	 */
	bool operator!=(const double & real) const;	
	
	/**
	 * Méthode pour écrire un nombre complex dans un flux sortant.
	 *
	 * @param out le flux dans lequel on souhaite écrire le nombre complexe.
	 *
	 * @see operator<<( ostream & out , const complex &).
	 */
	void write(ostream & out) const;

};
	
/**
 * Opérateur pour l'écriture d'un complexe dans un flux.
 *
 * @param out le flux dans lequel on souhaite écrire le nombre complexe.
 *
 * @param c le complex que l'on souohaite écrire dans le flux.
 *
 * @return une référence vers le flux dans lequel le complexe a été écrit.
 */
extern ostream & operator<<(ostream & out , const complex & c);

/**
 * Opérateur pour faire l'addition d'un réel et d'un complexe
 *
 * @param real le nombre réel.
 *
 * @param c un complexe.
 *
 * @return le résultat de l'addition du complexe et du réel.
 */
complex operator+(const double & real , const complex & c);

/**
 * Opérateur pour faire soustraction d'un réel et d'un complexe
 *
 * @param real le nombre réel.
 *
 * @param c un complexe.
 *
 * @return le résultat de la soustraction de real et de c.
 */
complex operator-(const double & real , const complex & c);

/**
 * Opérateur pour faire division d'un réel par un complexe
 *
 * @param real le nombre réel.
 *
 * @param c un complexe.
 *
 * @return le résultat de la division de real par c.
 */
complex operator/(const double & real , const complex & c);

/**
 * Opérateur pour faire la multiplication d'un réel et d'un complexe
 *
 * @param real le nombre réel.
 *
 * @param c un complexe.
 *
 * @return le résultat de la multiplication du complexe et du réel.
 */
complex operator*(const double & real , const complex & c);

/**
 * Opérateur unaire -
 *
 * @param c un complexe dont on veut prendre l'opposé.
 *
 * @return l'opposé du nombre complexe passé en argument.
 */
complex operator - (const complex & c);


#endif
