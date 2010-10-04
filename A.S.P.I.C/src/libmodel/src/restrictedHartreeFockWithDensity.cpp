/* 
* The model library of the A.S.P.I.C. 
 * Written and directed by François Lodier suport.aspic@gmail.com.
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
#include "restrictedHartreeFockWithDensity.h"
#include <generalizedEigenSolver.h>

//////////////////////////////////////////////////////////////////////////////////////////
// constructor of the class.
//////////////////////////////////////////////////////////////////////////////////////////
restrictedHartreeFockWithDensity::restrictedHartreeFockWithDensity(void)
: restrictedHartreeFock() , Density2CoefficientsIsUp2Date(true)
{
	;
}

//////////////////////////////////////////////////////////////////////////////////////////
// destructor of the class.
//////////////////////////////////////////////////////////////////////////////////////////
restrictedHartreeFockWithDensity::~restrictedHartreeFockWithDensity(void)
{
  clear();
}

//////////////////////////////////////////////////////////////////////////////////////////
// clear method of the class.
//////////////////////////////////////////////////////////////////////////////////////////
void restrictedHartreeFockWithDensity::clear(void)
{
	molecularSystem::clear();
	Density.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////
// density2CoefficientsIsUp2Date ?
//////////////////////////////////////////////////////////////////////////////////////////
bool restrictedHartreeFockWithDensity::density2CoefficientsIsUp2Date(void) const
{
	return Density2CoefficientsIsUp2Date;
}

//////////////////////////////////////////////////////////////////////////////////////////
// the GET method for the coefficients matrix.
//////////////////////////////////////////////////////////////////////////////////////////
matrixFull restrictedHartreeFockWithDensity::getCoefficients(void) const
{
	matrixFull coefficients;
	int row , column , nbrOfElectronsdiv2 , nbrOfBasisFunctions;
	
	// We just check that the converter is
	// up 2 date. If it's not then we do 
	// the update.
	if(!density2CoefficientsIsUp2Date()) {
		((restrictedHartreeFockWithDensity &)(*this)).updateCoefficients();
	}
	
	// Now we store the eigen vectors in the matrix
	// that will be returned.
	nbrOfBasisFunctions = getNbrOfBasisFunctions();
	nbrOfElectronsdiv2 = getNbrOfElectrons()/2;
	
	coefficients.setMatrixSize(getNbrOfBasisFunctions(),getNbrOfElectrons());
	
	for(row = 0 ; row < nbrOfBasisFunctions ; row++) {
		for(column = 0 ; column < nbrOfElectronsdiv2 ; column++) {
			coefficients(row,column) = getEigenVector(column)[row];
		}	
	}
	
	return coefficients;
}

//////////////////////////////////////////////////////////////////////////////////////////
// the GET method for the eigenvalues.
//////////////////////////////////////////////////////////////////////////////////////////
const containor<double> & restrictedHartreeFockWithDensity::getEigenValues(void) const
{
	// We just check that the converter is
	// up 2 date. If it's not then we do 
	// the update.
	if(!density2CoefficientsIsUp2Date()) {
		((restrictedHartreeFockWithDensity &)(*this)).updateCoefficients();
	}
	
	return Density2Coefficients.getEigenValues();
}

//////////////////////////////////////////////////////////////////////////////////////////
// the GET method for the eigenvalues.
//////////////////////////////////////////////////////////////////////////////////////////
const double & restrictedHartreeFockWithDensity::getEigenValue(const int & item) const
{
	assert(item >= 0);
	assert(item < getNbrOfBasisFunctions());
	
	// We just check that the converter is
	// up 2 date. If it's not then we do 
	// the update.
	if(!density2CoefficientsIsUp2Date()) {
		((restrictedHartreeFockWithDensity &)(*this)).updateCoefficients();
	}
	
	return Density2Coefficients.getEigenValue(item);
}

//////////////////////////////////////////////////////////////////////////////////////////
// the GET method for the eigenvalues.
//////////////////////////////////////////////////////////////////////////////////////////
const vector & restrictedHartreeFockWithDensity::getEigenVector(const int & item) const
{
	assert(item >= 0);
	assert(item < getNbrOfBasisFunctions());
	
	// We just check that the converter is
	// up 2 date. If it's not then we do 
	// the update.
	if(!density2CoefficientsIsUp2Date()) {
		((restrictedHartreeFockWithDensity &)(*this)).updateCoefficients();
	}
	
	return Density2Coefficients.getEigenVector(item);
}

//////////////////////////////////////////////////////////////////////////////////////////
// the GET method for the eigenvalues.
//////////////////////////////////////////////////////////////////////////////////////////
const containor<vector> & restrictedHartreeFockWithDensity::getEigenVectors(void) const
{
	// We just check that the converter is
	// up 2 date. If it's not then we do 
	// the update.
	if(!density2CoefficientsIsUp2Date()) {
		((restrictedHartreeFockWithDensity &)(*this)).updateCoefficients();
	}
	
	return Density2Coefficients.getEigenVectors();
}

//////////////////////////////////////////////////////////////////////////////////////////
// the GET method for the density matrix.
//////////////////////////////////////////////////////////////////////////////////////////
const matrixSymetric & restrictedHartreeFockWithDensity::getDensity(void) const
{
	return Density;
}

//////////////////////////////////////////////////////////////////////////////////////////
// the SET method for the coefficients.
//////////////////////////////////////////////////////////////////////////////////////////
void restrictedHartreeFockWithDensity::setCoefficients(const matrixFull & coefficients)
{
	assert(coefficients.getNbrOfRows() == getNbrOfBasisFunctions());
	assert(coefficients.getNbrOfColumns() == (getNbrOfElectrons() /2) );

	int nbrOfBasisFunctions , nbrOfElectronsdiv2 , i , j , k;
	
	nbrOfBasisFunctions = getNbrOfBasisFunctions();
	nbrOfElectronsdiv2 = getNbrOfElectrons() /2;
	
	for(i=0 ; i < nbrOfBasisFunctions ; i++) {
		for(j=i ; j < nbrOfBasisFunctions ; j++) {
			Density(i,j) = 0;
			for(k=0 ; k < nbrOfElectronsdiv2 ; k++) {
				Density(i,j) += coefficients(i,k) * coefficients(j,k);
			}
		} 
	}		

	Density2CoefficientsIsUp2Date = false;
}

//////////////////////////////////////////////////////////////////////////////////////////
// the SET method for the coefficients.
//////////////////////////////////////////////////////////////////////////////////////////
void restrictedHartreeFockWithDensity::setDensity(const matrixSymetric & density)
{
	assert(density.getNbrOfRows() == getNbrOfBasisFunctions());
	Density = density;
	Density2CoefficientsIsUp2Date = false;
}

////////////////////////////////////////////////////////////////////////////////////////////
// This is the method SET molecularSystem that is reimplemented from here
////////////////////////////////////////////////////////////////////////////////////////////
void restrictedHartreeFockWithDensity::setMolecularSystem(const molecule & mol , bool showProgress)
{
	
	if(mol.getNbrOfElectrons() % 2) {
		cerr << "Error : in void restrictedHartreeFockWithDensity::setMolecularSystem(const molecule & mol , bool showProgress)" << endl;
		cerr << "Error : the number of electrons of the molecule is not odd." << endl;
		cerr << "Error : unable to build the molecular system" << endl;
		clear();
		return;
	}
	Density.setMatrixSize(mol.getNbrOfBasisFunctions());
	molecularSystem::setMolecularSystem(mol,showProgress);
	Density2CoefficientsIsUp2Date = false;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Method that updates the coefficients.
////////////////////////////////////////////////////////////////////////////////////////////
void restrictedHartreeFockWithDensity::updateCoefficients(void) 
{
	// If we are up to date, then we do nothing.
	if(density2CoefficientsIsUp2Date())
		return;

	Density2Coefficients.solve(getFockMatrix() ,getOverlap());
	Density2CoefficientsIsUp2Date = true;
}
