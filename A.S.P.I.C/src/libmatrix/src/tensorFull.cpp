/* 
* The martix library of the A.S.P.I.C. 
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
#include "tensorFull.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The default constructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tensorFull::tensorFull(void)
: tensorNonSparse()
{
	setTensorType(FULL);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The constructor with the sizes.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tensorFull::tensorFull(const int & nbrOfRows , const int & nbrOfColumns , const int & nbrOfWidths , const int & nbrOfHeights)
: tensorNonSparse()
{
	setTensorType(FULL);
	setTensorSize( nbrOfRows , nbrOfColumns , nbrOfWidths , nbrOfHeights);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The destructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tensorFull::~tensorFull(void)
{
	;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Method that computes the total storage size.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int tensorFull::computeStorageSize(void) const
{
	return getNbrOfRows() * getNbrOfColumns() * getNbrOfWidths() * getNbrOfHeights();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Method that computes the position in the storage array of the element row x column.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int tensorFull::hash(const int & row , const int & column , const int & width , const int & height) const
{
	assert( column >= 0 );
	assert( column < getNbrOfColumns() );

	assert( height >= 0 );
	assert( height < getNbrOfHeights() );

	assert( row >= 0 );
	assert( row < getNbrOfRows() );
		
	assert( width >= 0 );
	assert( width < getNbrOfWidths() );

	return height + ( width + (column + row*getNbrOfColumns()) * getNbrOfWidths()) * getNbrOfHeights();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that writes down the tensor in a human readable way
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void tensorFull::writePretty(ostream & outStream , int width , int precision) const
{
	int i , j , k , l , nbrOfColumns , nbrOfRows , nbrOfWidths , nbrOfHeights;	
	
	
	nbrOfRows = getNbrOfRows();
	nbrOfColumns = getNbrOfColumns();
	nbrOfWidths = getNbrOfWidths();
	nbrOfHeights =	getNbrOfHeights();
	
	for(i=0 ; i < nbrOfRows ; i++) {
		for(j=0 ; j < nbrOfColumns ; j++) {
			
		
			for(k=0 ; k < nbrOfWidths ; k++) {
				for(l=0 ; l < nbrOfHeights ; l++) {
					outStream << setw(width) << setprecision(precision) << getCoefficient(i,j , k , l);
				}
			}
			
			outStream << endl;
		}
	}
}
