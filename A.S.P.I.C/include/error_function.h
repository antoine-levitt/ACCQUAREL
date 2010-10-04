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
#ifndef _ERROR_FUNCTION_
#define _ERROR_FUNCTION_

/**
 * @file error_function.h
 * Ceci permet de définir la fonction d'erreur.
 * En effet les plateformes windows ne possèdent pas de fonction
 * d'erreur dans leur librairies standard cette implémentation a 
 * donc été portée afon que le module GimC fonctionne sous windows.
 */

/**
	 * Méthode pour calculer la valeur de la fonction d'erreur
	 * \f[
	 * erf(x) = \frac{2}{\sqrt{\pi}} \int_{0}^{x}{\exp(-u^2)du}
	 * \f].
	 *
	 * La fonction d'erreur est de la forme 
	 *
	 * @param x la borne supérieure de l'intégrale.
	 * @return la valeur de la fonction d'erreur au point \f$x\f$.
	 */
	extern double error_function(double x);
	

	/**
	 * Méthode pour calculer la valeur de la fonction d'erreur complémentaire
	 * \f[
	 * erfc(x) = \frac{2}{\sqrt{\pi}} \int_{x}^{+\infty}{\exp(-u^2)du}
	 * \f].
	 *
	 * @param x la borne supérieure de l'intégrale.
	 * @return la valeur de la fonction d'erreur au point.
	 */
	extern double error_function_c(double x);

#endif


