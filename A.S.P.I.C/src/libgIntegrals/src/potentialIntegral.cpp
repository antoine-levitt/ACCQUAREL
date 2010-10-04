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
#include "potentialIntegral.h"

#include <assert.h>
#include "hermiteProjector.h"
#include "hermiteIntegrationBuffer.h"
#include <ipoint.h>


//////////////////////////////////////////
// Constructeur par défaut.
//////////////////////////////////////////
potentialIntegral::potentialIntegral(void)
{
	;
}

//////////////////////////////////////////////////////////////////
// Constructeur de copie interdit.
//////////////////////////////////////////////////////////////////
//potentialIntegral::potentialIntegral(const potentialIntegral & qe)
//{
//	assert(0);
//}

///////////////////////////////////////////////////////////////////////////////////////
// Méthode qui calcule l'intégrale de potential pour deux polynomes gaussiens.
////////////////////////////////////////////////////////////////////////////////////////
double potentialIntegral::computePotential(const gaussianPolynome3D & gp_a , const gaussianPolynome3D & gp_b , const dpoint<3> & center)
{
	// Les projections des polynomes gaussiens
	// sur la base de Hermite qui convient.
	hermiteProjector h_ab;
	
	// Le buffer avec les valeurs d'intégrales.
	hermiteIntegrationBuffer buffer;
	
	// Des varialbles pour mettre des résultats
	// intermédiaires.
	double value=0 , coefficient;
	
	// On construit le produit des polynomes
	// sur la base des polynomes de Hermite.
	h_ab.multiply(gp_a,gp_b);
	
	// On construit le buffer.
	buffer.setCenter(computePotentialCenter(h_ab,center));
	buffer.setExponent(computePotentialExponent(h_ab,center));
	buffer.setIntegrationBufferSize(h_ab.getHermiteDegreeMax());
	
	// On assemble le résultat du buffer pour les 
	// projections des polynomes.
	for(ipoint<3> i = h_ab.getHermiteFirst() ; i != h_ab.getHermiteEnd() ; i = h_ab.getHermiteNext(i)) {
			value += h_ab.getHermiteCoeffcient(i) * buffer.getIntegrationValue(i);
	}	
	
	// On caclule le coefficient devant :
	coefficient = computePotentialCoefficient(h_ab , center);
	
	// On renvoie la valeur.
	return coefficient* h_ab.getCoefficient() * value;
}


//////////////////////////////////////////////////////////////////////////////////////
// Méthode qui calcule le centrage pour l'intégrale de potential.
//////////////////////////////////////////////////////////////////////////////////////
dpoint<3> potentialIntegral::computePotentialCenter(const hermiteProjector & h_ab , const dpoint<3> & center) 
{
	dpoint<3> potentialCenter;
	potentialCenter = center - h_ab.getCenter();

	if(potentialCenter.norme_2() < 1E-10) {
		potentialCenter.setDatas(0);
	}
	
	return potentialCenter;
}

////////////////////////////////////////////////////////////////////////////////////
// Lorsque nous calculons une intégrales de potential, il sort un coefficient de la 
// forme 2 * Pi ^(5/2) * alpha * beta * sqrt( alpha + beta).   
// En fait cette méthode calcul ce coefficient à partir des données alpha et beta.
/////////////////////////////////////////////////////////////////////////////////////
double potentialIntegral::computePotentialCoefficient(const hermiteProjector & h_ab , const dpoint<3> & center)
{
	double coefficient , alpha;
	
	//On récupère les informations.
	alpha = h_ab.getExponent();
	
	// Calcul du coefficient.
	coefficient =  2 * M_PI / alpha;
	
	// on renvoie ce qu'il faut.
	return coefficient;
}

///////////////////////////////////////////////////////////////////////////////////////
// Méthode qui calcul ce qui sert d'exposant dans l'intégrale
// de potential.
////////////////////////////////////////////////////////////////////////////////////////
double potentialIntegral::computePotentialExponent(const hermiteProjector & h_ab , const dpoint<3> & center)
{
	double exponent , alpha;
	
	// On récupère ce qu'il faut.
	alpha = h_ab.getExponent();

	// On calcule la valeur.
	exponent = alpha;

	// On renvoie la valeur.
	return exponent;
}



///////////////////////////////////////////////////////////////////////////////
// Méthode qui calcule l'intégrale pour deux fonctions de base et un centre.
////////////////////////////////////////////////////////////////////////////////
double potentialIntegral::getPotentialValue(const gaussianBasisFunction & Phi_A , const gaussianBasisFunction & Phi_B , dpoint<3> center, const atom::distanceUnit & unit) const
{
	int i , j  , Ka , Kb ;
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


	///////////////////////////////////////////////////////////////////////////
	// Gestion de l'unité dexpression des coordonées du centre :
	//
	//	
	///////////////////////////////////////////////////////////////////////////
	if(unit != atom::ATOMIC_UNIT) {
		
		if(unit == atom::ANGSTROEM) {
			center = atom::angstroem2AtomicUnitConverter(center);	
		} else {
			assert(0);
		}

	}
	
	// On boucle sur le nombre de contractions de la premiÃ¨re
	// fonction de base.
	for(i=0 ; i < Ka ; i++) {
		gp_a = Phi_A.getGaussianPolynome3D(i);

		// On boucle ensuite sur le nombre de recouvrement
		// de la seconde fonction de base.
		for(j=0 ; j < Kb ; j++) {
			gp_b = Phi_B.getGaussianPolynome3D(j);
			value += computePotential(gp_a,gp_b,center);
		}
	}

	// Il ne reste plus qu'a retourner la valeure.
	return value;

}

///////////////////////////////////////////////////////////////////////////////
// Méthode qui calcule l'intégrale pour les derivees de deux fonctions de base et un centre.
////////////////////////////////////////////////////////////////////////////////
double potentialIntegral::getPotentialValue(const gaussianBasisFunction & Phi_A , const ipoint<3> & deg_a, const gaussianBasisFunction & Phi_B , const ipoint<3> & deg_b, dpoint<3> center, const atom::distanceUnit & unit) const
{
	int Ka , Kb ;
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

	Ka = Phi_A.getNbrOfContractions();
	Kb = Phi_B.getNbrOfContractions();

	///////////////////////////////////////////////////////////////////////////
	// Gestion de l'unité dexpression des coordonées du centre :
	///////////////////////////////////////////////////////////////////////////
	if (unit != atom::ATOMIC_UNIT) {
	  if (unit == atom::ANGSTROEM) {
		center = atom::angstroem2AtomicUnitConverter(center);	
	  } else {
		assert(0);
	  }
	}
	
	// On boucle sur le nombre de contractions de la premiÃ¨re
	// fonction de base.
	for(int i=0 ; i < Ka ; i++) {
	  gp_a = Phi_A.getGaussianPolynome3D(i).derivate(deg_a);

	  // On boucle ensuite sur le nombre de recouvrement
	  // de la seconde fonction de base.
	  for(int j=0 ; j < Kb ; j++) {
		gp_b = Phi_B.getGaussianPolynome3D(j).derivate(deg_b);
		value += computePotential(gp_a,gp_b,center);
	  }
	}

	// Il ne reste plus qu'a retourner la valeur.
	return value;

}

