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
#include "hermiteIntegrationBuffer.h"
#define _USE_MATH_DEFINES
#include <math.h>

//////////////////////////////////////////////////////////////////////////
// Constructeur.
//////////////////////////////////////////////////////////////////////////
hermiteIntegrationBuffer::hermiteIntegrationBuffer(void)
{
	clearMonomeIntegrationBuffer();
}

//////////////////////////////////////////////////////////////////////////
// Méthode pour mettre le buffer à zéro.
//////////////////////////////////////////////////////////////////////////
void hermiteIntegrationBuffer::clearMonomeIntegrationBuffer(void)
{
	// On fait juste un clear du buffer.
	IntegrationBuffer.clear();
	
	//On dit que le buffer est à jour puisqu'il n'y a 
	// rien à calculer.
	IntegrationBufferIsUp2Date = true;
}

//////////////////////////////////////////////////////////////////////////
// Méthode pour calculer les intégrales (cas général).
//////////////////////////////////////////////////////////////////////////
void hermiteIntegrationBuffer::computeHermiteIntegrationBuffer(void)
{	
	// Lorsque le buffer est à jour je ne recalcule rien.
	if(IntegrationBufferIsUp2Date == true) {
		return;
	}
	
	// Lorsque nous avons un DegreeMax pourri il
	// n'y a rien à faire 
	if(getIntegrationBufferSize(0) == 0 || getIntegrationBufferSize(1) == 0 || getIntegrationBufferSize(2) == 0)
	{
		nothing2DoException();
	}
	
	////////////////////////////////////// La boucle pourrie ! //////////////////
	tpoint<hermiteBasisPolynome,3> hermiteBasisPolynomes;
	gaussianPolynome gPolynome;
	polynome px , pxy , pxyz;
	int rx , ry ,rz;

	// Initialisaiton du polynome gaussien.
	for(int i=0 ; i < 3 ; i++) {
		hermiteBasisPolynomes[i].setExponent(getExponent());
		hermiteBasisPolynomes[i].setCenter(getCenter()[i]);
	}

	// faut mettre le bon exposant pour la gaussienne.
	// sauf quand celui ci est nul, alors dans ce cas
	// on va juste rien faire du tout puisue l'on utilisera pas la gaussienne.
	if(getCenter() != dpoint<3>(0)) {
		gPolynome.setExponent(getExponent() * getCenter().sq_norme_2());
	}

	// On initialise le polynome de Hermite en Rz
	hermiteBasisPolynomes[0].setDegree(0);
	for( rx = 0 ;rx < getIntegrationBufferSize(0); rx++){


		// Remarque : Lorsque nous avons une composante 
		// du centre qui est nulle et un degré impaire 
		// pour celle ci alors la valeur est forcément 0.
		if( getCenter()[0] == 0 && (rx%2)) {
			for( ry=0 ; ry < getIntegrationBufferSize(1); ry++) {
				for( rz = 0 ; rz < getIntegrationBufferSize(2) ; rz++) {
					IntegrationBuffer(rx,ry,rz) = 0;
				}
			}

			// Si besoin est on construit le polynome
			// de hermite suivant.
			if(rx < (getIntegrationBufferSize(0)-1)) 
				hermiteBasisPolynomes[0].addDegree();

			// On skip la fin de la boucle.
			continue;
		}

		// on réinitialise le polynome de Hermite en Ry.
		hermiteBasisPolynomes[1].setDegree(0);
		for( ry=0 ; ry < getIntegrationBufferSize(1); ry++) {

			// Remarque : Lorsque nous avons une composante 
			// du centre qui est nulle et un degré impaire 
			// pour celle ci alors la valeur est forcément 0.
			if( getCenter()[1] == 0 && (ry%2)) {
				for( rz = 0 ; rz < getIntegrationBufferSize(2) ; rz++) {
					IntegrationBuffer(rx,ry,rz) = 0;
				}

				// Si besoin est on construit le polynome
				// de hermite suivant.
				if(ry < (getIntegrationBufferSize(1)-1)) 
					hermiteBasisPolynomes[1].addDegree();

				// On skip la fin de la boucle.	
				continue;
			}

			// Ici on calcule le produit intéermediare entre les polynomes avec Rx et Ry
			// pour ne pas avoir à la refaire tout le temps.
			pxy = hermiteBasisPolynomes[0].getPolynome() * hermiteBasisPolynomes[1].getPolynome();

			// On réinitialise le polynome de hermite en Rz.
			hermiteBasisPolynomes[2].setDegree(0);
			for( rz = 0 ; rz < getIntegrationBufferSize(2) ; rz++) {


				// Remarque : Lorsque nous avons une composante 
				// du centre qui est nulle et un degré impaire 
				// pour celle ci alors la valeur est forcément 0.
				if( getCenter()[2] == 0 && (rz%2)) {
					IntegrationBuffer(rx,ry,rz) = 0;
					
					// Si besoin est on construit le polynome
					// de hermite suivant.
					if(rz < (getIntegrationBufferSize(2)-1)) 
						hermiteBasisPolynomes[2].addDegree();

					// On skip la fin de la boucle.
					continue;
				}

				pxyz = pxy * hermiteBasisPolynomes[2].getPolynome();

				if(getCenter() != dpoint<3>(0)) {
					// Lorsque nous avons un centre non nul,
					// alors nous avons bien une intégrale de 
					// polynome gaussien à calculer.
					
					// On construit le polynome gaussien 
					// qui se trouve dans l'intégralle.
					gPolynome.setPolynomeCoefficients(pxyz.monomeMultiply(rx+ry+rz));

					// Ensuite on demande la primitive du polynome
					// gaussien ainsi construit.
					IntegrationBuffer(rx,ry,rz) =  gPolynome.primitive(1);
				} else {
					// sinon nous avons vraiment un truc tout con 
					// comme l'integrale de u^r entre 0 et 1 ...
					IntegrationBuffer(rx,ry,rz) =  pxyz.getData(0) / (rx+ry+rz+1.);
				}

				// Si besoin est on construit le polynome
				// de hermite suivant.
				if(rz < (getIntegrationBufferSize(2)-1)) 
					hermiteBasisPolynomes[2].addDegree();
			}

			// Si besoin est on construit le polynome
			// de hermite suivant.
			if(ry < (getIntegrationBufferSize(1)-1)) 
				hermiteBasisPolynomes[1].addDegree();
		}

		// Si besoin est on construit le polynome
		// de hermite suivant.
		if(rx < (getIntegrationBufferSize(0)-1)) 
			hermiteBasisPolynomes[0].addDegree();
	}

	// Faut dire que nous sommes à jour.
	IntegrationBufferIsUp2Date = true;
}

//////////////////////////////////////////////////////////////////////////
// Méthode pour connaitre le centre.
//////////////////////////////////////////////////////////////////////////
const dpoint<3> & hermiteIntegrationBuffer::getCenter(void) const
{
	return Center;
}

///////////////////////////////////////////////////////////////////////////
// Méthode pour connaitre le degré maximal.
///////////////////////////////////////////////////////////////////////////
const ipoint<3> hermiteIntegrationBuffer::getIntegrationBufferSize(void) const
{
	return IntegrationBuffer.getSizes();
}

///////////////////////////////////////////////////////////////////////////
// Méthode pour connaitre le degré maximal.
///////////////////////////////////////////////////////////////////////////
int hermiteIntegrationBuffer::getIntegrationBufferSize(const int & dimension) const
{
	assert(dimension >= 0);
	assert(dimension < 3);
	return IntegrationBuffer.getSize(dimension);
}

///////////////////////////////////////////////////////////////////////////
// Méthode pour connaitre le degré maximal.
///////////////////////////////////////////////////////////////////////////
double hermiteIntegrationBuffer::getExponent(void) const
{
	return Exponent;
}

////////////////////////////////////////////////////////////////////////////
// Méthode pour accéder aux valeurs de l'intégrale.
////////////////////////////////////////////////////////////////////////////
const double & hermiteIntegrationBuffer::getIntegrationValue(const int & rx , const int & ry , const int & rz) const
{
	// 1 vérification sur les degrés.
	assert(rx >= 0);
	assert(rx < getIntegrationBufferSize(0));
	assert(ry >= 0);
	assert(ry < getIntegrationBufferSize(1));
	assert(rz >= 0);
	assert(rz < getIntegrationBufferSize(2));

	if(IntegrationBufferIsUp2Date == false) {
		((hermiteIntegrationBuffer&)(*this)).computeHermiteIntegrationBuffer();
	}
	
	return IntegrationBuffer.getData(rx,ry,rz);
}

const double & hermiteIntegrationBuffer::getIntegrationValue(const ipoint<3> & r) const
{
	// 1 vérification sur les degrés.
	assert(r[0] >= 0);
	assert(r[0] < getIntegrationBufferSize(0));
	assert(r[1] >= 0);
	assert(r[1] < getIntegrationBufferSize(1));
	assert(r[2] >= 0);
	assert(r[2] < getIntegrationBufferSize(2));
	
	if(IntegrationBufferIsUp2Date == false) {
		((hermiteIntegrationBuffer&)(*this)).computeHermiteIntegrationBuffer();
	}
	
	// 2 on acède simplment au containor 3D.
	return IntegrationBuffer.getData(r[0],r[1],r[2]);
}

///////////////////////////////////////////////////////////////////////////
// Gestion de l'expetion ou il n'y a rien à faire.
///////////////////////////////////////////////////////////////////////////
void hermiteIntegrationBuffer::nothing2DoException(void) 
{
	cerr << "Warning : in void hermiteIntegrationBuffer::nothing2DoException(void)" << endl;
	cerr << "Warning : no data will be computed." << endl;
	IntegrationBufferIsUp2Date = true;
}

//////////////////////////////////////////////////////////////////////////
// Méthode pour calculer le point de centrage 
// de nos intégrales.
///////////////////////////////////////////////////////////////////////////
void hermiteIntegrationBuffer::setCenter(const dpoint<3> & center)
{
	Center =  center;
	IntegrationBufferIsUp2Date = false;
}

////////////////////////////////////////////////////////////////////////////////
// M2thode pour calculer l'expsoant.
////////////////////////////////////////////////////////////////////////////////
void hermiteIntegrationBuffer::setExponent(double exposant)
{
	assert(exposant > 0);
	Exponent =  exposant;
	IntegrationBufferIsUp2Date = false;
}

////////////////////////////////////////////////////////////////////////////////
// Méthode pour les taille.
////////////////////////////////////////////////////////////////////////////////
void hermiteIntegrationBuffer::setIntegrationBufferSize(ipoint<3> degree)
{
	for(int i=0 ; i < 3 ; i++) {
		// On vérifie que le degré est positif.
		assert(degree[i] >=0);
		degree[i]++;
	}
	
	IntegrationBuffer.setSizes(degree);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// La classe pour faire le calcul des polynome de hermite.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
//Le constructeur de la classe.
/////////////////////////////////////////////////////////////
hermiteBasisPolynome::hermiteBasisPolynome(void)
: Center(0) , 
Degree(-1) ,
Exponent(1) ,
HermiteBasisPolynomeIsInitialized(false) ,
HermiteBasisPolynomeIsUp2Date(false)
{
	;
}

////////////////////////////////////////////////////////////
// Destructeur.
////////////////////////////////////////////////////////////
hermiteBasisPolynome::~hermiteBasisPolynome(void)
{
	;
}

////////////////////////////////////////////////////////////
// méthode pour augmenter le degré.
////////////////////////////////////////////////////////////
void hermiteBasisPolynome::addDegree(void)
{
	setDegree(getDegree()+1);	
}


////////////////////////////////////////////////////////////
// Méthode qui calcule le polynome de Hermite souhaité.
////////////////////////////////////////////////////////////
void hermiteBasisPolynome::computeHermiteBasisPolynome(void)
{	
	// Lorsque nous sommes déja à jour, alors nous n'avons rien
	// à faire. On sort d'ici tout de suite.
	if(HermiteBasisPolynomeIsUp2Date) {
		return;
	}

	// Si le polynome n'est pas intialisé, ou si nous voulons 
	// un degré qui est inférieur aue degré du polynome actuel
	// alors nous devons réinitialisé tout ...
	if( HermiteBasisPolynomeIsInitialized == false || getDegree() < HermiteBasisPolynome.getPolynomeDegree()) {
		initHermiteBasisPolynome();
	}

	//Dans le cas ou nous ne sommes pas avec un centre nul
	if(getCenter() != 0) {
		//Ensuite il ne reste plus qu'a augmenté le degré betement.
		while(HermiteBasisPolynome.getPolynomeDegree() < getDegree()) {
			HermiteBasisPolynome = HermiteBasisPolynome.derivate();
			((polynome &)HermiteBasisPolynome)/=(-getCenter());		
		}
	} else {
		computeHermiteBasisPolynomeNullCenterException(getDegree());
	}

	// On dit ensuite que nous sommes à jour.
	HermiteBasisPolynomeIsUp2Date = true;
}

///////////////////////////////////////////////////////////////
// Méthode qui calcule la tête du polynome dans le cas ou nous
// avons un centre nul.
///////////////////////////////////////////////////////////////
double __factoriel(int n , double result=1)
{
	assert( (n%2) == 0);
	if(n==0) {
		return result;
	} else {
		return __factoriel(n-2, result * 2 *(n-1) ) ;
	}
}

void hermiteBasisPolynome::computeHermiteBasisPolynomeNullCenterException(const int & degree)
{	
	// On vérifie aussi que l'initialisation à été 
	// correctement effectuée.
	assert(HermiteBasisPolynomeIsInitialized);
	
	// On vérifie simplement que si nous passons par là
	// que l'exception existe.
	assert(getCenter() == 0);


	// Nous sommes dans le cas particulier où le polynome
	// de Hermite n'est pas un polynome mais juste une constante.
	HermiteBasisPolynome.clearPolynome();	

	// Lorsque le degré est paire la valeur du polynome 
	// constant est v^r * H_r(0). Soit encore
	// (-2*exposant)^(r/2) * (r/2) !
	if((degree % 2)==0) {
		HermiteBasisPolynome.setPolynomeCoefficients(pow(-getExponent() , degree /2) * __factoriel(degree)) ;
	} 

	// lorsque le degré est impaire la valeur et zéro donc nous n'avons rien à faire.
}

////////////////////////////////////////////////////////////
// Méthode pour acéder au centre.
///////////////////////////////////////////////////////////
const double & hermiteBasisPolynome::getCenter(void) const
{
	return Center;
}

/////////////////////////////////////////////////////////////
// Méthode pour accéder au degré du polynome
/////////////////////////////////////////////////////////////
const int & hermiteBasisPolynome::getDegree(void) const
{
	return Degree;
}

////////////////////////////////////////////////////////////
// Méthode pour accéder à l'exposant.
////////////////////////////////////////////////////////////
const double & hermiteBasisPolynome::getExponent(void) const
{ 
	return Exponent;
}

/////////////////////////////////////////////////////////////
// Méthode pur accéder au polynome.
/////////////////////////////////////////////////////////////
const polynome & hermiteBasisPolynome::getPolynome(void) const
{
	if(HermiteBasisPolynomeIsUp2Date == false) {
		((hermiteBasisPolynome &)(*this)).computeHermiteBasisPolynome();
	}
	
	return HermiteBasisPolynome.getPolynomeCoefficients();
}

////////////////////////////////////////////////////////////////
// Méthode pour accéder aux coefficient du polynome.
////////////////////////////////////////////////////////////////
double hermiteBasisPolynome::getPolynomeCoefficient(const int & degree) const
{
	if(HermiteBasisPolynomeIsUp2Date == false) {
		((hermiteBasisPolynome &)(*this)).computeHermiteBasisPolynome();
	}
	
	return  HermiteBasisPolynome.getPolynomeCoefficient(degree);
}
//////////////////////////////////////////////////////////////////
// Méthode d'initialisation de la gaussienne.
//////////////////////////////////////////////////////////////////
void hermiteBasisPolynome::initHermiteBasisPolynome(void)
{
	
	// On met ce qu'il faut comme Centre :
	if(getCenter() != 0) {
		HermiteBasisPolynome.setExponent( getExponent() * pow(getCenter(),2));
	} else {
		// normalement on s'en fout !!
		// puisque l'on utilisera les polynomes de Hermites
		// que lorsque l'on en aura besoin !
		;
	}
	
	// on réinitialise le polynome à 1 et on dit que le degré est zéro.
	HermiteBasisPolynome.setPolynomeCoefficients(1);
	Degree = 0;

	// On dit que l'initialisation à été faite.
	HermiteBasisPolynomeIsInitialized =true;
}

//////////////////////////////////////////////////////////////////
// Méthode pour fixer le centre
//////////////////////////////////////////////////////////////////
void hermiteBasisPolynome::setCenter(const double & center)
{
	Center = center;
	HermiteBasisPolynomeIsInitialized = false;
	HermiteBasisPolynomeIsUp2Date = false;
}

//////////////////////////////////////////////////////////////////
// Méthode pour fixer un degré.
//////////////////////////////////////////////////////////////////
void hermiteBasisPolynome::setDegree(const int & degree)
{
	Degree = degree;
	HermiteBasisPolynomeIsUp2Date = false;
}

////////////////////////////////////////////////////////////////
// Méthode pour fixer l'exposant.
/////////////////////////////////////////////////////////////////
void hermiteBasisPolynome::setExponent(const double & exponent)
{
	assert(exponent > 0);
	Exponent = exponent;
	HermiteBasisPolynomeIsInitialized = false;
	HermiteBasisPolynomeIsUp2Date = false;
}

/////////////////////////////////////////////////////////////////
// Méthode pour écrire le résultat dans un flux.
//////////////////////////////////////////////////////////////////
void hermiteBasisPolynome::write(ostream & out) const
{
	out << getPolynome();
}

//////////////////////////////////////////////////////////////////////
// Opérateur externe pour écire un hermiteBasisPolynome dans un flux.
//////////////////////////////////////////////////////////////////////
ostream & operator<<(ostream & out , const hermiteBasisPolynome & hBP)
{
	hBP.write(out);
	return out;
}

