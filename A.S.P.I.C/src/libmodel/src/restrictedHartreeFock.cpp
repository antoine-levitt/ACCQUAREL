/* 
 * The model library of A.S.P.I.C. 
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
#include "restrictedHartreeFock.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor of the class. 
//////////////////////////////////////////////////////////////////////////////////////////////////
restrictedHartreeFock::restrictedHartreeFock(void)
: molecularSystem()
{
	;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Destructor of the class. 
//////////////////////////////////////////////////////////////////////////////////////////////////
restrictedHartreeFock::~restrictedHartreeFock(void)
{
	;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// This method computes the RHF Energy.
//////////////////////////////////////////////////////////////////////////////////////////////////
matrixSymetric restrictedHartreeFock::applyTensor4RHF(const matrixSymetric & density) const 
{
	assert(density.getNbrOfRows() == getNbrOfBasisFunctions());
	
	matrixSymetric fockMatrix;
	fockMatrix = 2 * applyCoulomb(density);
	fockMatrix -= applyExchange(density);
	return fockMatrix;
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// This method computes the RHF Energy.
//////////////////////////////////////////////////////////////////////////////////////////////////
double restrictedHartreeFock::getElectronsEnergy(const matrixSymetric & density) const
{
	
	assert(density.getNbrOfRows() == getNbrOfBasisFunctions());
	
	int i,k , nbrOfRows;
	double energy;
	matrixSymetric energyMatrix;
	
	// le 2 a cause du spin: chaque couche est supposee contenir deux electrons,
	// un de spin up et l'autre de spin down (c'est la signification du mot Restricted dans RHF)
	// au lieu de mettre h::D, on va calculer 2 h::D
	energyMatrix = 2*getHamilton();
	// a cause du spin encore, il ne faut pas calculer 1/2 D::A::D, mais autre chose
	// sans spin, A = Coulomb - Exchange
	// avec spin, A = 4 Coulomb - 2 Exchange, donc il faut calculer
	// 1/2 D::(4 Coulomb - 2 Exchange)::D = 2 D::Coulomb::D - D::Exchange::D
	//                                    = (2 D::Coulomb - D::Exchange)::D
	// applyCoulomb renvoie D::Coulomb, et applyExchange renvoie D::Exchange
	energyMatrix += applyTensor4RHF(density);

	// Pour toutes ces histoires de Spin, cf. these E. Cances, p37-39.
	
	/*energyMatrix += 2*applyCoulomb(density);
	energyMatrix -= applyExchange(density);*/

	
	nbrOfRows = density.getNbrOfRows();
	energy = 0;
	for(i = 0 ; i < nbrOfRows ; i++) {
		for(k=0 ; k < nbrOfRows ; k++) {
			energy += energyMatrix(i,k) * density(k,i);
		}
	}

	return energy;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// This method computes the RHF Energy.
//////////////////////////////////////////////////////////////////////////////////////////////////
double restrictedHartreeFock::getElectronsEnergy(void) const
{
	return getElectronsEnergy(getDensity());
}	

//////////////////////////////////////////////////////////////////////////////////////////////////
// This method computes the RHF Energy.
//////////////////////////////////////////////////////////////////////////////////////////////////
double restrictedHartreeFock::getEnergy(const matrixSymetric & density) const
{
	return getElectronsEnergy(density) + getNucleiEnergy();
}	


//////////////////////////////////////////////////////////////////////////////////////////////////
// This method computes the RHF Energy.
//////////////////////////////////////////////////////////////////////////////////////////////////
double restrictedHartreeFock::getEnergy(void) const
{
	return getEnergy(getDensity());
}	

//////////////////////////////////////////////////////////////////////////////////////////////////
// Method that builds the fock matrix for an Hartree Fock model. 
//////////////////////////////////////////////////////////////////////////////////////////////////
matrixSymetric restrictedHartreeFock::getFockMatrix(const matrixSymetric & density) const
{
	matrixSymetric fockMatrix;

	// fockMatrix = hamilton + 2 * J(D) - K(D).
	// d'après le handbook donc on applique cela betement. 
	fockMatrix = getHamilton();
	//fockMatrix += 2 * applyCoulomb(density);
	//fockMatrix -= applyExchange(density);
	fockMatrix += applyTensor4RHF(density);
	return fockMatrix;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Method that builds the fock matrix for an Hartree Fock model. 
//////////////////////////////////////////////////////////////////////////////////////////////////
matrixSymetric restrictedHartreeFock::getFockMatrix(void) const
{
	return getFockMatrix(getDensity());
}
