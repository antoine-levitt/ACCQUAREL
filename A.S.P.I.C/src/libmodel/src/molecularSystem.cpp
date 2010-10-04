/* 
* The model library of the A.S.P.I.C. 
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
#include "molecularSystem.h"
#include <moleculeMapper.h>
#include <coulombIntegral.h>
#include <potentialIntegral.h>
#include <monoElectronicIntegral.h>
#include <time.h>

///////////////////////////////////////////////////////////////////////////////////
// Constructor of the class molecule.
///////////////////////////////////////////////////////////////////////////////////
molecularSystem::molecularSystem(void)
:Coulomb() , 
CoulombComputationTime(0) ,  
Hamilton() , 
HamiltonComputationTime(0),
NbrOfElectrons(0) ,
NucleiEnergy(0) ,
Overlap() ,
OverlapComputationTime(0)
{
	;
}

///////////////////////////////////////////////////////////////////////////////////
// Destructor of the class molecule.
///////////////////////////////////////////////////////////////////////////////////
molecularSystem::~molecularSystem(void)
{
	;
}

////////////////////////////////////////////////////////////////////////////////////
// Method that applies the coulomb tensor to a density matrix.
////////////////////////////////////////////////////////////////////////////////////
matrixSymetric molecularSystem::applyCoulomb(const matrixSymetric & density) const
{

	assert(density.getNbrOfRows() == getNbrOfBasisFunctions());
	
	matrixSymetric coulombMatrix;
	
	if(empty()) {
		return coulombMatrix;
	} 
	
	int i , j , k , l , nbrOfRows = getNbrOfBasisFunctions();
	double value;
	
	coulombMatrix.setMatrixSize(getNbrOfBasisFunctions());
	
	for(i=0 ; i < nbrOfRows ; i++) {
		for(j=i ; j < nbrOfRows ; j++) {
		
			value = 0;
			
			for(k=0 ; k < nbrOfRows ; k++) {
				value += Coulomb.getCoefficient(i,j,k,k) * density(k,k);
				for(l=k+1 ; l < nbrOfRows ; l++) {
					value += 2*Coulomb.getCoefficient(i,j,k,l) * density.getCoefficient(k,l);
				}		
			}
			
			coulombMatrix.setCoefficient(i,j,value);
		}		
	}
	
	return coulombMatrix;
}

////////////////////////////////////////////////////////////////////////////////////
// Method that applies the exchange tensor to a density matrix.
////////////////////////////////////////////////////////////////////////////////////
matrixSymetric molecularSystem::applyExchange(const matrixSymetric & density) const
{
	
	assert(density.getNbrOfRows() == getNbrOfBasisFunctions());
	
	matrixSymetric exchangeMatrix;
	
	if(empty()) {
		return exchangeMatrix;
	}
	
	int i , j , k , l , nbrOfRows = getNbrOfBasisFunctions();
	double value;
	
	exchangeMatrix.setMatrixSize(nbrOfRows);
	
	for(i=0 ; i < nbrOfRows ; i++) {
		for(j=i ; j < nbrOfRows ; j++) {
			
			value = 0;
			for(k=0 ; k < nbrOfRows ; k++) {
				for(l=0 ; l < nbrOfRows ; l++) {
					value += Coulomb.getCoefficient(i,k,j,l) * density.getCoefficient(k,l);
				}		
			}
			
			exchangeMatrix.setCoefficient(i,j,value);
		}		
	}
	
	return exchangeMatrix;
}


///////////////////////////////////////////////////////////////////////////////////
// Method clear for the molecular system class
///////////////////////////////////////////////////////////////////////////////////
void molecularSystem::clear(void)
{
	Coulomb.clear();
	Hamilton.clear();
	Overlap.clear();
	NbrOfElectrons = 0;
	NucleiEnergy = 0;
}


////////////////////////////////////////////////////////////////////////////////////
// Method that computes (fills) the coulomb tensor.
////////////////////////////////////////////////////////////////////////////////////
void molecularSystem::computeCoulomb(const molecule & mol , bool showProgress)
{
	assert(getNbrOfBasisFunctions() == mol.getNbrOfBasisFunctions());
	
	
	int  i, j , k , l , nbrOfBasisFunctions;
	moleculeMapper map;
	coulombIntegral coulomb;
	gaussianBasisFunction phi_a , phi_b , phi_c , phi_d;
	double value=0 , done , nextMark=0;
	time_t start,end;

	// From now on we have started the coulomb computation.
	time(&start);

	map.attach(mol);
	nbrOfBasisFunctions = getNbrOfBasisFunctions();
	
	
	if(showProgress) {
		clog << "Log : molecularSystem::computeCoulomb   |*" << ends;
		nextMark = 10;
	}

	for(i=0 ; i < nbrOfBasisFunctions ; i++) {
			phi_a = map.getBasisFunction4Global(i);				
			
			for(j=i ; j < nbrOfBasisFunctions ; j++) {
				phi_b = map.getBasisFunction4Global(j);

				// this is to see the progression of the computation.
				if(showProgress) {
					done =  100. * double(i*nbrOfBasisFunctions + j ) / double(nbrOfBasisFunctions * nbrOfBasisFunctions);
					if(done >= nextMark) {
						clog << "-" << ends;
						nextMark += 10;
					}
				}

				for(k=i ; k < nbrOfBasisFunctions ; k++) {
					phi_c = map.getBasisFunction4Global(k);
					
					for(l=k ; l < nbrOfBasisFunctions ; l++) {
						
						
						if(k==i && l<j) {
							continue;
						} 
						
						phi_d = map.getBasisFunction4Global(l);
						value = coulomb.getCoulombValue(phi_a, phi_b , phi_c,phi_d);
						
						Coulomb.setCoefficient(i,j,k,l,value);
					} // fin du for l.
				} // fin du for k.	
			} // fin du for j.
		} // fin du for i.
	

	// This is now over...
	// we store the enlapsed time in the coulombComputationTime
	// member.

	time(&end);
	CoulombComputationTime = difftime(end,start);
	
	if(showProgress) {
		cout << "-*| done in " << CoulombComputationTime << " second(s)" << endl;
	}

}

////////////////////////////////////////////////////////////////////////////////////
// Method that computes (fills) the hamilton matrix.
////////////////////////////////////////////////////////////////////////////////////
void molecularSystem::computeHamilton(const molecule & mol , bool showProgress)
{

	assert(getNbrOfBasisFunctions() == mol.getNbrOfBasisFunctions());

	int nbrOfBasisFunctions , i, j;
	moleculeMapper map;
	potentialIntegral potential;
	monoElectronicIntegral overlap;
	gaussianBasisFunction phi_a , phi_b;
	double value , done , nextMark=0;
	time_t start,end;

	// From now on we have started the coulomb computation.
	time(&start);

	if(showProgress) {
		clog << "Log : molecularSystem::computeHamilton  |*" << ends;
		nextMark = 10;
	}

	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();
	
		
	
	for(i=0 ; i < nbrOfBasisFunctions ; i++) {		
		phi_a = map.getBasisFunction4Global(i);
		
		for(j=i ; j < nbrOfBasisFunctions ; j++) {
			
			phi_b = map.getBasisFunction4Global(j);
			
			value = overlap.getKineticValue(phi_a,phi_b) / 2.;
			
			for(int atom = 0 ; atom < mol.getNbrOfAtoms() ; atom++) {
				value -= mol.getNbrOfProtons4Atom(atom) * potential.getPotentialValue(phi_a,phi_b,mol.getPosition4Atom(atom) , mol.getDistanceUnit4Atom(atom));					 
			}
			
			Hamilton.setCoefficient(i,j,value);

				// this is to see the progression of the computation.
				if(showProgress) {
					done =  100. * double(i*nbrOfBasisFunctions + j ) / double(nbrOfBasisFunctions * nbrOfBasisFunctions);
					if(done >= nextMark) {
						clog << '-' << ends;
						nextMark += 10;
					}
				}

		}	
	}

	// This is now over...
	// we store the enlapsed time in the coulombComputationTime
	// member.
	time(&end);
	HamiltonComputationTime = difftime(end,start);

	if(showProgress) {
		cout << "-*| done in " << HamiltonComputationTime << " second(s)" << endl;
	}
}	

////////////////////////////////////////////////////////////////////////////////////
// Method that the nuclei energy.
////////////////////////////////////////////////////////////////////////////////////
void molecularSystem::computeNucleiEnergy(const molecule & mol)
{
	int k , l , nbrOfAtoms;
	dpoint<3> xk , xl;
	
	NucleiEnergy=0;
	nbrOfAtoms = mol.getNbrOfAtoms();
	
	for(k=0 ; k < nbrOfAtoms ; k++) {
	for(l=k+1 ; l < nbrOfAtoms ; l++) {
	xk = mol.getAtom(k).getPosition(atom::ATOMIC_UNIT);
	xl = mol.getAtom(l).getPosition(atom::ATOMIC_UNIT);			
	NucleiEnergy += mol.getAtom(k).getNbrOfProtons() * mol.getAtom(l).getNbrOfProtons() / (xk-xl).norme_2();
	}		
}
	
}

////////////////////////////////////////////////////////////////////////////////////
// Method that computes (fills) the overlap matrix.
////////////////////////////////////////////////////////////////////////////////////
void molecularSystem::computeOverlap(const molecule & mol , bool showProgress)
{
	assert(getNbrOfBasisFunctions() == mol.getNbrOfBasisFunctions());
	
	int nbrOfBasisFunctions , i, j;
	moleculeMapper map;
	monoElectronicIntegral overlap;
	gaussianBasisFunction phi_a , phi_b;
	double value, done , nextMark=0;
	time_t start,end;

	// From now on we have started the coulomb computation.
	time(&start);
	
	if(showProgress) {
		clog << "Log : molecularSystem::computeOverlap   |*" << ends;
		nextMark = 10;
	}
	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();
		
	for(i=0 ; i < nbrOfBasisFunctions ; i++) {		
		phi_a = map.getBasisFunction4Global(i);
		
		for(j=i ; j < nbrOfBasisFunctions ; j++) {
			phi_b = map.getBasisFunction4Global(j);
			
			value = overlap.getOverlapValue(phi_a,phi_b);
						
			Overlap.setCoefficient(i,j,value);

			// this is to see the progression of the computation.
				if(showProgress) {
					done =  100. * double(i*nbrOfBasisFunctions + j ) / double(nbrOfBasisFunctions * nbrOfBasisFunctions);
					if(done >= nextMark) {
						clog << '-' << ends;
						nextMark += 10;
					}
				}

		}	
	}
	
	// This is now over...
	// we store the enlapsed time in the coulombComputationTime
	// member.
	time(&end);
	OverlapComputationTime =  difftime(end,start);

		if(showProgress) {
		cout << "-*| done in " << OverlapComputationTime << " second(s)" << endl;
	}

}

///////////////////////////////////////////////////////////////////////////////////
// Method empty for the molecular system class
///////////////////////////////////////////////////////////////////////////////////
bool molecularSystem::empty(void) const
{
	if(Coulomb.empty()) {
		return true;
	} else {
		return false;
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the coulomb tensor.
///////////////////////////////////////////////////////////////////////////////////////
const tensorSymetric & molecularSystem::getCoulomb(void) const
{
	return Coulomb;
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the coulomb computaiton time.
///////////////////////////////////////////////////////////////////////////////////////
const double & molecularSystem::getCoulombComputationTime(void) const
{
	return CoulombComputationTime;
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the Hamiltonian matrix.
///////////////////////////////////////////////////////////////////////////////////////
const matrixSymetric & molecularSystem::getHamilton(void) const
{
	return Hamilton;
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the hamilton computation time.
///////////////////////////////////////////////////////////////////////////////////////
const double & molecularSystem::getHamiltonComputationTime(void) const
{
	return HamiltonComputationTime;
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the hamilton computation time.
///////////////////////////////////////////////////////////////////////////////////////
const double & molecularSystem::getNucleiEnergy(void) const
{
	return NucleiEnergy;
}


///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the Overlap matrix.
///////////////////////////////////////////////////////////////////////////////////////
const matrixSymetric & molecularSystem::getOverlap(void) const
{
	return Overlap;
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the overlap computation time.
///////////////////////////////////////////////////////////////////////////////////////
const double & molecularSystem::getOverlapComputationTime(void) const
{
	return OverlapComputationTime;
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the  number of basis functions.
///////////////////////////////////////////////////////////////////////////////////////
const int & molecularSystem::getNbrOfBasisFunctions(void) const
{
	return Coulomb.getNbrOfRows();
}

///////////////////////////////////////////////////////////////////////////////////////
// Method GET for the  number of electrons.
///////////////////////////////////////////////////////////////////////////////////////
const int & molecularSystem::getNbrOfElectrons(void) const
{
	return NbrOfElectrons;
}

////////////////////////////////////////////////////////////////////////////////////////
// Method SET from a molecule.
////////////////////////////////////////////////////////////////////////////////////////
void molecularSystem::setMolecularSystem(const molecule & mol , bool showProgress)
{
	setNbrOfBasisFunctions(mol.getNbrOfBasisFunctions());
	setNbrOfElectrons(mol.getNbrOfElectrons());
	computeCoulomb(mol , showProgress);
	computeHamilton(mol , showProgress);
	computeNucleiEnergy(mol);
	computeOverlap(mol , showProgress);
}

////////////////////////////////////////////////////////////////////////////////////////
// Method SET for the number of basis functions.
////////////////////////////////////////////////////////////////////////////////////////
void molecularSystem::setNbrOfBasisFunctions(const int & nbrOfBasisFunctions)
{
	assert(nbrOfBasisFunctions > 0);
	
	Coulomb.setTensorSize(nbrOfBasisFunctions);
	Hamilton.setMatrixSize(nbrOfBasisFunctions);
	Overlap.setMatrixSize(nbrOfBasisFunctions);
}

////////////////////////////////////////////////////////////////////////////////////////
// Method SET for the number of electrons.
////////////////////////////////////////////////////////////////////////////////////////
void molecularSystem::setNbrOfElectrons(const int & nbrOfElectrons)
{
	assert(nbrOfElectrons >= 0);
	
	NbrOfElectrons = nbrOfElectrons;
}