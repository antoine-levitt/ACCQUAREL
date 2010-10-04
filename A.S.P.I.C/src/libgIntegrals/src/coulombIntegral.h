/* 
 * The gaussian integrals library of A.S.P.I.C. 
 * Written and directed by François Lodier support.aspic@gmail.com.
 * Class modified by Frederic Legoll Nov 06
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
#ifndef _COULOMB_INTEGRAL_
#define _COULOMB_INTEGRAL_

#include <dpoint.h>
#include <gaussianBasisFunction.h>
#include "hermiteProjector.h"

/**
 * Classe pour calculer les intégrales de Coulomb.
 */
class coulombIntegral 
{
 private:
  
 protected:
/**
 * Méthode pour calculer le point de centrage des intégrales de couclomb.
 */
  static dpoint<3> computeCoulombCenter(const hermiteProjector & h_ab , const hermiteProjector & h_cd);
	
/**
 * Méthode pour calculer le coefficient à mettre devant les intégrales de coulomb.
 */
  static double computeCoulombCoefficient(const hermiteProjector & h_ab , const hermiteProjector & h_cd);
		
/**
 * Méthode pour calculer l'exposant de l'intégrale de coulomb.
 */
  static double computeCoulombExponent(const hermiteProjector & h_ab , const hermiteProjector & h_cd);
		
/**
 * Méthode pour calculer l'intégrale bielectronique pour quatre polynomes gaussiens.
 */
  static double computeCoulomb(const gaussianPolynome3D & gp_a , const gaussianPolynome3D & gp_b , const gaussianPolynome3D & gp_c , const gaussianPolynome3D & gp_d);
		
 public:
/**
 * Constructeur par défaut.
 */
  coulombIntegral(void);
		
  /**
   * Méthode pour faire le calcul de l'intégrale bi électronique pour quatre fonctions de base.
   */
  double getCoulombValue(const gaussianBasisFunction & phi_a , const gaussianBasisFunction & phi_b , const gaussianBasisFunction & phi_c , const gaussianBasisFunction & phi_d) const;

  /**
   * Méthode pour faire le calcul de l'intégrale bi électronique pour les derivees de quatre fonctions de base.
   */
  double getCoulombValue(const gaussianBasisFunction & phi_a , const ipoint<3> & deg_a , const gaussianBasisFunction & phi_b , const ipoint<3> & deg_b , const gaussianBasisFunction & phi_c , const ipoint<3> & deg_c, const gaussianBasisFunction & phi_d, const ipoint<3> & deg_d) const;

  
	
};

#endif


