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
#ifndef _HERMITE_PROJECTOR_
#define _HERMITE_PROJECTOR_


#include <containor.h>
#include <gaussian3D.h>
#include <ipoint.h>
#include <polynome3D.h>
#include <tpoint.h>

/**
 * Classe des polynomes de Hermite.
 */
class hermiteProjector :  public gaussian3D
{
	private:	
		
		/**
		 * 
		 */
		polynome3D HermitePolynome;

		/**
		 * This buffer contains the monomes on the hermite 
		 * Basis.
		 */
		tpoint<containor<polynome3D> ,3> MonomeBuffer;
		
		/**
		 * Boolean value to know if computation of the buffer 
		 * needs to be performed.
		 */
		bool MonomeBufferIsComputed;

		/**
		 * Boolean to know if the initialisation of the buffer
		 * needs to be performed.
		 */
		bool MonomeBufferIsInitialized;

		/**
		 * The maximal size for the monome buffer.
		 */
		ipoint<3> MonomeBufferSize; 
	
	protected:
	
		/**
		 * Method GET for the Monome decomposition on the
		 * Hermite polynome basis.
		 */
		polynome3D getMonomeDecomposition(const ipoint<3> & degree) const;
		
		/**
		 * Method to initialize the monome decomposition.
		 */
		void initMonomeBuffer(void);

		/**
		 * Method that mutliplys by a single varible an Hermite polynome.
		 */
		polynome3D shiftHermitePolynome(const polynome3D & hermitePolynome , const int & dim) const;

		/**
		 * Method that clears the monome buffer.
		 */
		void clearMonomeBuffer(void);

		/**
		 * Method that computes the decompositions for all the mmonomes.
		 */
		void computeMonomeDecomposition(void);
	
		/**
		 * Method to find the decompostion of x^degree
		 * on the hermite polynome basis.
		 */
		const polynome3D & getMonomeBufferData(const int & dimension ,const int & degree) const;


		/**
		 * Method get to acces the greatest degree of the buffer.
		 */
		const ipoint<3> & getMonomeBufferSize(void) const;

		/**
		 * Method get to acces the greatest degree of the buffer.
		 *
		 * If dimension == 0 returns max degree for x.
		 * If dimension == 1 returns max degree for y.
		 * If dimension == 2 returns max degree for z.
		 */
		const int & getMonomeBufferSize(const int & dimension) const;
		
		/**
		 * Method SET for the monome buffer size.
		 */
		void setMonomeBufferSize(ipoint<3> degreeMax);

		/**
		 * Method SET for the monome buffer size.
		 */
		void setMonomeBufferSize(const int & deg_x , const int & deg_y , const int & deg_z);

	public:

		/**
		 * Constructor.
		 */
		hermiteProjector(void);

		/**
		 * Method GET to access the coefficients of the Hermite polynome.
		 */
		double getHermiteCoeffcient(const ipoint<3> & degree) const;

		/**
		 * Methode GET for the end of the Hermite Polynome.
		 */
		ipoint<3> getHermiteEnd(void) const;
		
		/**
		 * Method GET for the begining of the Hermite Polynome.
		 */
		ipoint<3> getHermiteFirst(void) const;

		/**
		 * Method GET to know the degrees following the monome
		 * with degree degree.
		 */
		ipoint<3> getHermiteNext(const ipoint<3> & degree) const;
		
		/**
		 * Method get for the max degree of the polynome.
		 */
		ipoint<3> getHermiteDegreeMax(void) const;

		/**
		 * Method that multiplies two gaussian polynomes.
		 * The result is stored as an Hermite polynome in 
		 * the calling object.
		 */
		void multiply(const gaussianPolynome3D & gp_1 , const gaussianPolynome3D & gp_2);

		/**
		 * Method that converts a gaussian polynome
		 * onto an Hermite polynome.
		 */
		void set(const gaussianPolynome3D & gp);

		/**
		 * Method SET for the gaussian exponent.
		 *
		 * @warning the value of the exponent shall be 
		 * a stricly positive one.
		 */
		void setExponent(const double & exponent);
};

#endif

