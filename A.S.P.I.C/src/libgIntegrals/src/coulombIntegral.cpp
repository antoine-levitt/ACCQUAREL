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
#include "coulombIntegral.h"

#include <assert.h>
#include "hermiteProjector.h"
#include "hermiteIntegrationBuffer.h"
#include <ipoint.h>


coulombIntegral::coulombIntegral(void)
{
	;
}

//coulombIntegral::coulombIntegral(const coulombIntegral & qe)
//{
//	assert(0);
//}

///////////////////////////////////////////////////////////////////////////////////////
// Méthode qui calcule l'intégrale de coulomb pour deux polynomes gaussiens.
////////////////////////////////////////////////////////////////////////////////////////
double coulombIntegral::computeCoulomb(const gaussianPolynome3D & gp_a , const gaussianPolynome3D & gp_b , const gaussianPolynome3D & gp_c , const gaussianPolynome3D & gp_d)
{
	// Les projections des polynomes gaussiens
	// sur la base de Hermite qui convient.
	hermiteProjector h_ab;
	hermiteProjector h_cd;
	
	// Le buffer avec les valeurs d'intégrales.
	hermiteIntegrationBuffer buffer;
	
	// Des varialbles pour mettre des résultats
	// intermédiaires.
	double value=0 , signe , coefficient;
	
	// On construit le produit des polynomes
	// sur la base des polynomes de Hermite.
	h_ab.multiply(gp_a,gp_b);
	h_cd.multiply(gp_c,gp_d);
	
	// On construit le buffer.
	buffer.setCenter(computeCoulombCenter(h_ab,h_cd));
	buffer.setExponent(computeCoulombExponent(h_ab,h_cd));
	buffer.setIntegrationBufferSize(h_ab.getHermiteDegreeMax() + h_cd.getHermiteDegreeMax());
	
	// On assemble le résultat du buffer pour les 
	// projections des polynomes.
	for(ipoint<3> i = h_ab.getHermiteFirst() ; i != h_ab.getHermiteEnd() ; i = h_ab.getHermiteNext(i)) {
		for(ipoint<3> j = h_cd.getHermiteFirst() ; j != h_cd.getHermiteEnd() ; j = h_cd.getHermiteNext(j)) {
			if(j.norme_1() % 2) {
				signe = -1;
			}	else {
				signe = 1;
			}
			value += signe * h_ab.getHermiteCoeffcient(i) * h_cd.getHermiteCoeffcient(j) * buffer.getIntegrationValue(i+j);
		}
	}	
	
	// On caclule le coefficient devant :
	coefficient = computeCoulombCoefficient(h_ab , h_cd);
	
	// On renvoie la valeur.
	return coefficient* h_ab.getCoefficient() * h_cd.getCoefficient() * value;
}


//////////////////////////////////////////////////////////////////////////////////////
// Méthode qui calcule le centrage pour l'intégrale de coulomb.
//////////////////////////////////////////////////////////////////////////////////////
dpoint<3> coulombIntegral::computeCoulombCenter(const hermiteProjector & h_ab , const hermiteProjector & h_cd) 
{
	dpoint<3> coulombCenter;
	coulombCenter = h_cd.getCenter() - h_ab.getCenter();

	if(coulombCenter.norme_2() < 10E-10) {
		coulombCenter.setDatas(0);
	}

	return coulombCenter;
}

////////////////////////////////////////////////////////////////////////////////////
// Lorsque nous calculons une intégrales de coulomb, il sort un coefficient de la 
// forme 2 * Pi ^(5/2) * alpha * beta * sqrt( alpha + beta).   
// En fait cette méthode calcul ce coefficient à partir des données alpha et beta.
/////////////////////////////////////////////////////////////////////////////////////
double coulombIntegral::computeCoulombCoefficient(const hermiteProjector & h_ab , const hermiteProjector & h_cd)
{
	double coefficient , alpha , beta;
	
	//On récupère les informations.
	alpha = h_ab.getExponent();
	beta  = h_cd.getExponent();
	
	
	// Calcul du coefficient.
	coefficient =  2 * pow(M_PI , 5 / 2.) / (alpha * beta * sqrt(alpha +beta) );
	
	// on renvoie ce qu'il faut.
	return coefficient;
}

///////////////////////////////////////////////////////////////////////////////////////
// Méthode qui calcul ce qui sert d'exposant dans l'intégrale
// de coulomb.
////////////////////////////////////////////////////////////////////////////////////////
double coulombIntegral::computeCoulombExponent(const hermiteProjector & h_ab , const hermiteProjector & h_cd)
{
	double exponent , alpha , beta;
	
	// On récupère ce qu'il faut.
	alpha = h_ab.getExponent();
	beta  = h_cd.getExponent();

	// On calcule la valeur.
	exponent = alpha * beta / (alpha + beta);

	// On renvoie la valeur.
	return exponent;
}


///////////////////////////////////////////////////////////////////////////////////////
// Méthode pour le calcul de l'integrale bielectronique avec 4 fonctions de base
// 
////////////////////////////////////////////////////////////////////////////////////////
double coulombIntegral::getCoulombValue(const gaussianBasisFunction & Phi_A, const gaussianBasisFunction & Phi_B, const gaussianBasisFunction & Phi_C, const gaussianBasisFunction & Phi_D) const
{
	int i , j , k , l , Ka , Kb , Kc, Kd;
	gaussianPolynome3D gp_a, gp_b , gp_c , gp_d;
	double value;

	/////////////////////////////////////////////////
	// - Initialisation de la valeur de l'intÃ©gralle
	// Ã  zÃ©ro.
	// - Initialisation des nombres de contractions afin 
	// de ne pas avoir Ã  rappeler tout le temps
	// la mÃ©thode get_nbr_of_contraction.
	/////////////////////////////////////////////////
	value = 0;

	Ka = Phi_A.getNbrOfContractions();
	Kb = Phi_B.getNbrOfContractions();
	Kc = Phi_C.getNbrOfContractions();
	Kd = Phi_D.getNbrOfContractions();


	// On boucle sur le nombre de contractions de la premiÃ¨re
	// fonction de base.
	for(i=0 ; i < Ka ; i++) {
		gp_a = Phi_A.getGaussianPolynome3D(i);

		// On boucle ensuite sur le nombre de recouvrement
		// de la seconde fonction de base.
		for(j=0 ; j < Kb ; j++) {
			gp_b = Phi_B.getGaussianPolynome3D(j);

			for(k=0 ; k < Kc ; k++) {
				gp_c = Phi_C.getGaussianPolynome3D(k);
				
				for(l=0 ; l < Kd ; l++) {
					gp_d = Phi_D.getGaussianPolynome3D(l);
					value += computeCoulomb(gp_a,gp_b,gp_c,gp_d);
				}
			}
		}
	}

	// Il ne reste plus qu'a retourner la valeure.
	return value;

}

///////////////////////////////////////////////////////////////////////////////////////
// Méthode pour le calcul de l'integrale bielectronique pour les derivees de quatre fonctions de base.
// 
////////////////////////////////////////////////////////////////////////////////////////
double coulombIntegral::getCoulombValue(const gaussianBasisFunction & phi_a , const ipoint<3> & deg_a , const gaussianBasisFunction & phi_b , const ipoint<3> & deg_b , const gaussianBasisFunction & phi_c , const ipoint<3> & deg_c, const gaussianBasisFunction & phi_d, const ipoint<3> & deg_d) const
{
  int Ka , Kb , Kc, Kd;
  gaussianPolynome3D gp_a, gp_b , gp_c , gp_d;
  double value;

  /////////////////////////////////////////////////
  // - Initialisation de la valeur de l'intÃ©gralle
  // Ã  zÃ©ro.
  // - Initialisation des nombres de contractions afin 
  // de ne pas avoir Ã  rappeler tout le temps
  // la mÃ©thode get_nbr_of_contraction.
  /////////////////////////////////////////////////
  value = 0.;

  Ka = phi_a.getNbrOfContractions();
  Kb = phi_b.getNbrOfContractions();
  Kc = phi_c.getNbrOfContractions();
  Kd = phi_d.getNbrOfContractions();


  // On boucle sur le nombre de contractions de la premiÃ¨re
  // fonction de base.
  for(int i=0 ; i < Ka ; i++) {
	gp_a = phi_a.getGaussianPolynome3D(i).derivate(deg_a);

	// On boucle ensuite sur le nombre de recouvrement
	// de la seconde fonction de base.
	for(int j=0 ; j < Kb ; j++) {
	  gp_b = phi_b.getGaussianPolynome3D(j).derivate(deg_b);

	  for(int k=0 ; k < Kc ; k++) {
		gp_c = phi_c.getGaussianPolynome3D(k).derivate(deg_c);
				
		for(int l=0 ; l < Kd ; l++) {
		  gp_d = phi_d.getGaussianPolynome3D(l).derivate(deg_d);
		  value += computeCoulomb(gp_a,gp_b,gp_c,gp_d);
		}
	  }
	}
  }

  // Il ne reste plus qu'a retourner la valeur.
  return value;
}
