/* 
* The gaussian integrals library of A.S.P.I.C. 
 * Written and directed by Fran�ois Lodier support.aspic@gmail.com.
 *
 * Copyright (C) 2005  Fran�ois Lodier
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
#include "hermiteProjector.h"

/////////////////////////////////////////////////////////////////////////////////////
// Constructeur.
/////////////////////////////////////////////////////////////////////////////////////
hermiteProjector::hermiteProjector(void)
:  
MonomeBufferIsComputed(false),
MonomeBufferIsInitialized(false),
MonomeBufferSize(0)
{
	;
}

/////////////////////////////////////////////////////////////////////////////////////
// Méthode clear pour le buffer.
/////////////////////////////////////////////////////////////////////////////////////
void hermiteProjector::clearMonomeBuffer(void)
{
	for(int i=0 ; i<3 ; i++) {
		MonomeBuffer[i].clear();
		MonomeBufferSize[i] = 0;
	}
	MonomeBufferIsComputed = false;
	MonomeBufferIsInitialized = false;
}

/////////////////////////////////////////////////////////////////////////////////////
// Méthode pour calculer toutes les décompoistions de 
// monome sur la base des polynomes de Hermite.
/////////////////////////////////////////////////////////////////////////////////////
void hermiteProjector::computeMonomeDecomposition(void)
{	
	int dim , degree;

	// Lorsque le buffer est up 2 date alors on peut se permettre
	// de ne rien faire.
	if(MonomeBufferIsComputed) {
		return;
	}
	
	// Initialisation.
	initMonomeBuffer();

	// Récrence.
	for(dim = 0  ; dim < 3 ; dim++) {

		// En fait on utilise betement le fait
		// x^degree = x * x^(degree-1);
		//          = x * (un polynome de hermite connu). 
		for(degree = 1 ; degree < getMonomeBufferSize(dim) ; degree ++) {
				// On shift simplement le polynome d'avant.
				MonomeBuffer[dim][degree] = shiftHermitePolynome( MonomeBuffer[dim][degree-1] , dim);
		} // Fin du for degré.
	} // Fin du  for dim.
		
	// Il ne reste plus qu'à dire que nous avons fait notre calcul.
	MonomeBufferIsComputed = true;	
}

/////////////////////////////////////////////////////////////////////////////////////
// Méthode pour connaitre le degré du polynome de Hermite.
/////////////////////////////////////////////////////////////////////////////////////
ipoint<3> hermiteProjector::getHermiteDegreeMax(void) const
{
	return HermitePolynome.getUniformMax4Degree();
}

/////////////////////////////////////////////////////////////////////////////////////
// Méthode pour connaitre la taille du buffer.
/////////////////////////////////////////////////////////////////////////////////////
const ipoint<3> & hermiteProjector::getMonomeBufferSize(void) const
{
	return MonomeBufferSize;
}

///////////////////////////////////////////////////////////////////////////////////////
// M�thode pour connaitre la taille du buffer dans une direction 
// donnée.
///////////////////////////////////////////////////////////////////////////////////////
const int & hermiteProjector::getMonomeBufferSize(const int & dimension) const
{
	// Vérification d'usage.
	assert(dimension >= 0);
	assert(dimension < 3);

	// On retourne le degré max dans la direction
	return getMonomeBufferSize()[dimension];
}

/////////////////////////////////////////////////////////////////////////////
// Méthode pour acc�der aux coeffcient des polynomes de Hermite.
/////////////////////////////////////////////////////////////////////////////
double hermiteProjector::getHermiteCoeffcient(const ipoint<3> & degree) const
{
	return HermitePolynome.getData(degree);
}

//////////////////////////////////////////////////////////////////////////////
// Méthode pour connaitre le degré marqueur de fin du polynome de
// Hermite.
//////////////////////////////////////////////////////////////////////////////
ipoint<3> hermiteProjector::getHermiteEnd(void) const
{
	return HermitePolynome.end();
}

//////////////////////////////////////////////////////////////////////////////
// Méthode pour connaitre le premier degré du polynome de Hermite.
///////////////////////////////////////////////////////////////////////////////
ipoint<3> hermiteProjector::getHermiteFirst(void) const
{
	return HermitePolynome.begin();
}

////////////////////////////////////////////////////////////////////////////////
// Méthode pour connaitre le degré suivant du polynome de Hermite.
////////////////////////////////////////////////////////////////////////////////
ipoint<3> hermiteProjector::getHermiteNext(const ipoint<3> & degree) const
{
	return HermitePolynome.next(degree);
}

///////////////////////////////////////////////////////////////////////////////////////
// Méthode pour accéder à la décomposition d'un monome sur la base des polynomes
// de Hermite.
////////////////////////////////////////////////////////////////////////////////////////
const polynome3D & hermiteProjector::getMonomeBufferData(const int & dimension , const int & degree) const
{
	// 0 - Vérifications.
	assert(dimension >= 0);
	assert(dimension < 3);
	assert( degree < getMonomeBufferSize(dimension));

	// 1 - Bon Que doit on faire : il faut vérifier 
	// le Buffer est up2date ...
	if(MonomeBufferIsInitialized == false || MonomeBufferIsComputed == false) {
		cerr << "Error : in hermiteProjector::getMonome4Dimension( ... )" << endl;
		cerr << "Error : the monome decomposition buffer is marked as deprecated." << endl;
		cerr << "Error : result may vary." << endl;
	}
	
	return MonomeBuffer[dimension][degree];
}

///////////////////////////////////////////////////////////////////////////////////////
// Méthode pour accéder à la décomposition d'un monome sur la base des polynomes
// de Hermite.
////////////////////////////////////////////////////////////////////////////////////////
polynome3D hermiteProjector::getMonomeDecomposition(const ipoint<3> & degree) const
{
	// 0 - Vérifications.
	assert( degree[0] < getMonomeBufferSize(0));
	assert( degree[1] < getMonomeBufferSize(1));
	assert( degree[2] < getMonomeBufferSize(2));

	// 1 - On assemble le polynome pour
	// le degré demandé.
	polynome3D gp(1);
	for(int i=0 ; i < 3 ; i++) {
		gp *= getMonomeBufferData(i,degree[i]);
	}
	
	// 2 - On renvoie le polynome.
	return gp;
}

//////////////////////////////////////////////////////////////////////////////////////
// Méthode pour initiliser la décomposition des 
// monomes sur la base des polynomes de Hermite.
//////////////////////////////////////////////////////////////////////////////////////
void hermiteProjector::initMonomeBuffer(void)
{

	// Lorsque le membre MonomeDecompositionIsInitialized
	// est à true cela veut dire que nous n'avons rien à faire.
	if(MonomeBufferIsInitialized) {
		return;
	}

	// Sinon il faut initialiser notre récurence.
	// On utilise le Fait que sqrt(alpha)^0 *  H_0( sqrt(alpha)* u) = 1.
	int dim;
	ipoint<3> degree(0);

	for(dim = 0 ; dim < 3 ; dim++) {

		// On initialise 1 avec 1.
		MonomeBuffer[dim][0].setPolynomeCoefficients(1);
	
	} // fin du for dimension.

	// On dit maintenant que nous avons fait l'initialisation.
	MonomeBufferIsInitialized = true;
}

/////////////////////////////////////////////////////////////////////////////////////
// Multiplication de deux gaussiennes polynomes.
/////////////////////////////////////////////////////////////////////////////////////
void hermiteProjector::multiply(const gaussianPolynome3D & gp_a , const gaussianPolynome3D & gp_b)
{
	// 1. On commence par calculer le produit des deux
	// polynomes gaussiens de façon usuelle.
	gaussianPolynome3D gp = gp_a * gp_b;

	// 2. On projete ensuite le résultat sur l'objet 
	// courrant.
	set(gp);
}

/////////////////////////////////////////////////////////////////////////////////////
// M�thode set.
/////////////////////////////////////////////////////////////////////////////////////
void hermiteProjector::set(const gaussianPolynome3D & gp)
{
	
	// 1 - Construction de la mesure d'othogonalit�.
	// 1 - a . Sauvegarde du centre.
	// 1 - b . Sauvegarde du coefficient.
	// 1 - c . Mise à jour de l'exposant.
	setCenter(gp.getCenter());
	setCoefficient(gp.getCoefficient());
	setExponent(gp.getExponent());

	// 2 - Une fois nous connaissons la mesure
	// produit il faut projeter le polynome.
	
	// 2 - a . On commence par regarder un peu 
	// la taille des choses dont nous allons avoir
	// besoin.
	setMonomeBufferSize(gp.getPolynomeUniformMax4Degree());

	// 2 - b .  On calcule la projection de tout les
	// monomes.
	computeMonomeDecomposition();

	// 2 - c. On assemble ...
	ipoint<3> i;
	HermitePolynome.clear();
	for(i=gp.polynomeBegin() ; i != gp.polynomeEnd() ; i = gp.polynomeNext(i)) {
		HermitePolynome += gp.getPolynomeCoefficient(i) * getMonomeDecomposition(i);
	}

	// 3 - On nettoie le Buffer.
	clearMonomeBuffer();
}


/////////////////////////////////////////////////////////////////////////////////////
// Mise à la taille.
// Cette m�thode permet de sp�cifier les plus haut degr� pour 
// la d�composition de la base.
/////////////////////////////////////////////////////////////////////////////////////
void hermiteProjector::setMonomeBufferSize(ipoint<3> degreeMax)
{
	// On vérifie que le degré est bien tout positif.
	assert(degreeMax[0] >= 0);
	assert(degreeMax[1] >= 0);
	assert(degreeMax[2] >= 0);

	// On arrete tout de suite de parler en degr� mais
	// en nombre d'élément stock�s.
	int dim;
	for(dim = 0  ; dim < 3 ; dim++) {
		degreeMax[dim]++;
	}	

	// Lorsque nous sommes à la bonne taille on skip 
	// m�chament.
	if(degreeMax == getMonomeBufferSize()) {
		return;
	}	
	
	// Mise à la taille du buffer.
	for(dim = 0  ; dim < 3 ; dim++) {
		MonomeBuffer[dim].setSizes(degreeMax[dim] + 1);
	}	

	// On garde un trace des tailles.
	MonomeBufferSize = degreeMax;

	// Ici il vaut mieux réinitialiser les choses pour la récurence.
	MonomeBufferIsInitialized = false;
	
	// Ici il vaut mieux réinitialiser les choses pour la récurence.
	MonomeBufferIsComputed = false;
}

void hermiteProjector::setMonomeBufferSize(const int & deg_x, const int & deg_y , const int & deg_z)
{
	// On vérifie que le degré est bien tout positif.
	assert(deg_x >= 0);
	assert(deg_y >= 0);
	assert(deg_z >= 0);

	setMonomeBufferSize(ipoint<3>(deg_x,deg_y,deg_z));
}
///////////////////////////////////////////////////////////////////////////
// Méthode pour fixer la valeur de l'exposant de la gaussienne.
///////////////////////////////////////////////////////////////////////////
void hermiteProjector::setExponent(const double & exponent) 
{
	// Si on ne change pas la valeur on skippe.
	if(getExponent() == exponent) {
		return;
	}

	// Sinon, on copie la valeur
	gaussian3D::setExponent(exponent);

	// Il faut recalculer la récurence.
	MonomeBufferIsComputed = false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode de shiftage d'un polynome de Hermite.
//
// Description :
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
polynome3D hermiteProjector::shiftHermitePolynome(const polynome3D & hermitePolynome , const int & dim) const 
{
	
	// On vérifie que la dimension est la
	// dans la bonne plage de valeurs possibles.
	assert(dim >= 0);
	assert(dim < 3);

	ipoint<3> degree;
	double coefficient;
	polynome3D shiftPolynome;

	// On initialise le polynome shifté.
	shiftPolynome.clear();

	for(degree = hermitePolynome.begin() ; degree != hermitePolynome.end() ; degree = hermitePolynome.next(degree))
	{
		// On récupère le coefficient du polynome.
		coefficient = hermitePolynome.getData(degree);

		// La première partie en 
		// H_{degree + 1}(u) / (2 * alpha)
		degree[dim]++;
		shiftPolynome.setData(degree, shiftPolynome.getData(degree) + coefficient / (2 * getExponent()));
		degree[dim]--;

		// La première partie en 
		// 2 * degree * H_{degree - 1}(u)
		if(degree[dim] > 0) {
			degree[dim]--;
			shiftPolynome.setData(degree, shiftPolynome.getData(degree) +  (degree[dim] + 1) * coefficient );
			degree[dim]++;
		} 
	} // Fin du for degree.

	return shiftPolynome;
}


