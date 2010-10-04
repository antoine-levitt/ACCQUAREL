/* 
 * The gaussian integrals library of A.S.P.I.C. 
 * Written and directed by François Lodier support.aspic@gmail.com.
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
#include "monoElectronicIntegral.h"

#include <assert.h>

////////////////////////////////////////////////////////
// Constructeur de la classe.
////////////////////////////////////////////////////////
monoElectronicIntegral::monoElectronicIntegral(void)
{
	;
}

///////////////////////////////////////////////////////////////////////////////
// Fonctions pour le calcul d'intégrales intervenant dans le cas relativiste.
// (dernier ajout le 23/09/2008 par G. Legendre (CEREMADE - université de Paris-Dauphine))
///////////////////////////////////////////////////////////////////////////////
double monoElectronicIntegral::computeDeriv(const gaussianPolynome3D & gp_a , const gaussianPolynome3D & gp_b, const int & dimension)
{
        gaussianPolynome3D gp;

        gp = gp_a.derivate(dimension) * gp_b;

        return gp.integral();
}

double monoElectronicIntegral::getDerivValue(const gaussianBasisFunction & Phi_A , const gaussianBasisFunction & Phi_B, const int & dimension) const
{
        int i , j , Ka , Kb;
        gaussianPolynome3D gp_a,gp_b;
        double value;

        value = 0;
        Ka = Phi_A.getNbrOfContractions();
        Kb = Phi_B.getNbrOfContractions();

        for(i=0 ; i < Ka ; i++) {
                gp_a = Phi_A.getGaussianPolynome3D(i);
                for(j=0 ; j < Kb ; j++) {
                        gp_b = Phi_B.getGaussianPolynome3D(j);
                        value += computeDeriv(gp_a,gp_b,dimension);
                }
        }
        return value;
}

double monoElectronicIntegral::computeXDeriv(const gaussianPolynome3D & gp_a , const gaussianPolynome3D & gp_b, const int & dimension1, const int & dimension2)
{
        gaussianPolynome3D gp;

        gp = gp_a.derivate(dimension1) * gp_b;

        return gp.centerMonomeMultiply(0.,dimension2).integral();
}

double monoElectronicIntegral::getXDerivValue(const gaussianBasisFunction & Phi_A , const gaussianBasisFunction & Phi_B, const int & dimension1, const int & dimension2) const
{
        int i , j , Ka , Kb;
        gaussianPolynome3D gp_a,gp_b;
        double value;

        value = 0;
        Ka = Phi_A.getNbrOfContractions();
        Kb = Phi_B.getNbrOfContractions();

        for(i=0 ; i < Ka ; i++) {
                gp_a = Phi_A.getGaussianPolynome3D(i);
                for(j=0 ; j < Kb ; j++) {
                        gp_b = Phi_B.getGaussianPolynome3D(j);
                        value += computeXDeriv(gp_a,gp_b,dimension1,dimension2);
                }
        }
        return value;
}
//////////////////
// Fin des ajouts.
//////////////////

////////////////////////////////////////////////////////////////////////////////
// Energie Cinétique pour deux polynomes gaussiens.
////////////////////////////////////////////////////////////////////////////////
double monoElectronicIntegral::computeKinetic(const gaussianPolynome3D & gp_a , const gaussianPolynome3D & gp_b)

{
	gaussianPolynome3D gp;
	double value = 0;

	for(int dimension = 0 ; dimension < 3 ; dimension++) {
		gp = gp_a.derivate(dimension) * gp_b.derivate(dimension);
		value += gp.integral();
	}

	return value;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// La méthode qui assemble l'énergie cinétique pour deux fonctions
// de base.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double monoElectronicIntegral::getKineticValue(const gaussianBasisFunction & Phi_A , const gaussianBasisFunction & Phi_B) const
{
	int i , j , Ka , Kb;
	gaussianPolynome3D gp_a, gp_b;
	double value;

	/////////////////////////////////////////////////
	// - Initialisation de la valeur de l'intégralle
	// à zéro.
	// - Initialisation des nombres de contractions afin 
	// de ne pas avoir à rappeler tout le temps
	// la méthode get_nbr_of_contraction.
	/////////////////////////////////////////////////

	value = 0;
	Ka = Phi_A.getNbrOfContractions();
	Kb = Phi_B.getNbrOfContractions();

	/////////////////////////////////////////////////////
	// On boucle sur le nombre de contractions de la première
	// fonction de base.
	/////////////////////////////////////////////////////
	for(i=0 ; i < Ka ; i++) {
		gp_a = Phi_A.getGaussianPolynome3D(i);
		/////////////////////////////////////////////////////
		// On boucle ensuite sur le nombre de recouvrement 
		// de la seconde fonction de base.
		//////////////////////////////////////////////////////
		for(j=0 ; j < Kb ; j++) {
			gp_b = Phi_B.getGaussianPolynome3D(j);
			value += computeKinetic(gp_a,gp_b);	
		}
	}
	///////////////////////////////////////////////////
	// Il ne reste plus qu'a retourner la valeure.
	///////////////////////////////////////////////////
	return value;
}	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// La méthode qui assemble l'énergie cinétique pour la derivee de deux fonctions
// de base.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double monoElectronicIntegral::getKineticValue(const gaussianBasisFunction & Phi_A , const ipoint<3> & deg_a, const gaussianBasisFunction & Phi_B , const ipoint<3> & deg_b) const
{
	int Ka , Kb;
	gaussianPolynome3D gp_a, gp_b;
	double value;
	
	/////////////////////////////////////////////////
	// - Initialisation de la valeur de l'intégralle
	// à zéro.
	// - Initialisation des nombres de contractions afin 
	// de ne pas avoir à rappeler tout le temps
	// la méthode get_nbr_of_contraction.
	/////////////////////////////////////////////////
	
	value = 0.;
	Ka = Phi_A.getNbrOfContractions();
	Kb = Phi_B.getNbrOfContractions();
	
	/////////////////////////////////////////////////////
	// On boucle sur le nombre de contractions de la première
	// fonction de base.
	/////////////////////////////////////////////////////
	for(int i=0 ; i < Ka ; i++) {
		gp_a = Phi_A.getGaussianPolynome3D(i).derivate(deg_a);
		/////////////////////////////////////////////////////
		// On boucle ensuite sur le nombre de recouvrement 
		// de la seconde fonction de base.
		//////////////////////////////////////////////////////
		for(int j=0 ; j < Kb ; j++) {
			gp_b = Phi_B.getGaussianPolynome3D(j).derivate(deg_b);
			value += computeKinetic(gp_a,gp_b);	
		}
	}
	
	///////////////////////////////////////////////////
	// Il ne reste plus qu'a retourner la valeur.
	///////////////////////////////////////////////////
	return value;
}	

///////////////////////////////////////////////////////////
// La méthode qui calcule le recouvrement pour deux polynomes gaussiens.
///////////////////////////////////////////////////////////
double monoElectronicIntegral::computeOverlap(const gaussianPolynome3D & gp_a , const gaussianPolynome3D & gp_b)
{
	gaussianPolynome3D gp;	
	gp = gp_a * gp_b;
	return gp.integral();
}



///////////////////////////////////////////////////////////////
// La méthode qui assemble le recouvrement pour deux fonctions
// de base.
///////////////////////////////////////////////////////////////

double monoElectronicIntegral::getOverlapValue(const gaussianBasisFunction & Phi_A , const gaussianBasisFunction & Phi_B) const
{
	int i , j , Ka , Kb;
	gaussianPolynome3D gp_a, gp_b;
	double value;

	/////////////////////////////////////////////////
	// - Initialisation de la valeur de l'intégralle
	// à zéro.
	// - Initialisation des nombres de contractions afin 
	// de ne pas avoir à rappeler tout le temps
	// la méthode get_nbr_of_contraction.
	/////////////////////////////////////////////////
	value = 0;
	Ka = Phi_A.getNbrOfContractions();
	Kb = Phi_B.getNbrOfContractions();

	/////////////////////////////////////////////////////
	// On boucle sur le nombre de contractions de la première
	// fonction de base.
	/////////////////////////////////////////////////////
	for(i=0 ; i < Ka ; i++) {
		gp_a = Phi_A.getGaussianPolynome3D(i);

		/////////////////////////////////////////////////////
		// On boucle ensuite sur le nombre de recouvrement 
		// de la seconde fonction de base.
		//////////////////////////////////////////////////////
		for(j=0 ; j < Kb ; j++) {
			gp_b = Phi_B.getGaussianPolynome3D(j);
			value += computeOverlap(gp_a,gp_b);
		}
	}

	///////////////////////////////////////////////////
	// Il ne reste plus qu'a retourner la valeure.
	///////////////////////////////////////////////////
	return value;
}

///////////////////////////////////////////////////////////////
// La méthode qui assemble le recouvrement pour la derivee de deux fonctions
// de base.
///////////////////////////////////////////////////////////////
double monoElectronicIntegral::getOverlapValue(const gaussianBasisFunction & Phi_A , const ipoint<3> & deg_a , const gaussianBasisFunction & Phi_B , const ipoint<3> & deg_b) const
{
	int Ka , Kb;
	gaussianPolynome3D gp_a, gp_b;
	double value;
	
	/////////////////////////////////////////////////
	// - Initialisation de la valeur de l'intégralle
	// à zéro.
	// - Initialisation des nombres de contractions afin 
	// de ne pas avoir à rappeler tout le temps
	// la méthode get_nbr_of_contraction.
	/////////////////////////////////////////////////
	value = 0.;
	Ka = Phi_A.getNbrOfContractions();
	Kb = Phi_B.getNbrOfContractions();
	
	/////////////////////////////////////////////////////
	// On boucle sur le nombre de contractions de la première
	// fonction de base.
	/////////////////////////////////////////////////////
	for(int i=0 ; i < Ka ; i++) {
		gp_a = Phi_A.getGaussianPolynome3D(i).derivate(deg_a);
		
		/////////////////////////////////////////////////////
		// On boucle ensuite sur le nombre de recouvrement 
		// de la seconde fonction de base.
		//////////////////////////////////////////////////////
		for(int j=0 ; j < Kb ; j++) {
			gp_b = Phi_B.getGaussianPolynome3D(j).derivate(deg_b);
			value += computeOverlap(gp_a,gp_b);
		}
	}
	
	///////////////////////////////////////////////////
	// Il ne reste plus qu'a retourner la valeur.
	///////////////////////////////////////////////////
	return value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la normalisation d'une fonction de base.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void monoElectronicIntegral::normalize(gaussianBasisFunction & phi) const
{
	int i , K ;
	gaussianPolynome3D gp;
	double norm;

	K = phi.getNbrOfContractions();

	for(i=0 ; i < K ; i++) {
		gp = phi.getGaussianPolynome3D(i);
		gp.setCoefficient(1.);
		norm = sqrt(computeOverlap(gp,gp));
		phi.setCoefficient(i , phi.getCoefficient(i) / norm); 
	}	

}

