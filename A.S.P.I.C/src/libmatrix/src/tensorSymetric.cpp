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
#include "tensorSymetric.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The default constructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tensorSymetric::tensorSymetric(void)
: tensorNonSparse()
{
	setTensorType(SYMETRIC_FULL);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The constructor with the sizes.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tensorSymetric::tensorSymetric(const int & nbrOfRows)
: tensorNonSparse()
{
	setTensorType(SYMETRIC_FULL);
	setTensorSize(nbrOfRows);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The destructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tensorSymetric::~tensorSymetric(void)
{
	;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Method that computes the total storage size.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int tensorSymetric::computeStorageSize(void) const
{
	assert(getNbrOfColumns() == getNbrOfRows());
	assert(getNbrOfWidths() == getNbrOfHeights());
	assert(getNbrOfWidths() == getNbrOfRows());
	
	return ( getNbrOfRows() * (getNbrOfRows() + 1) / 2 ) * ( 1 + (getNbrOfRows() * (getNbrOfRows() + 1) / 2) ) / 2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Method that computes the position in the storage array of the element (row x column) x (width x height).
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int tensorSymetric::hash(const int & row , const int & column , const int & width , const int & height) const
{
	assert(getNbrOfRows() == getNbrOfColumns());
	assert(getNbrOfWidths() == getNbrOfHeights());
	assert(getNbrOfRows() == getNbrOfWidths());
		
	assert( column >= 0 );
	assert( column < getNbrOfColumns() );

	assert( height >= 0 );
	assert( height < getNbrOfHeights() );

	assert( row >= 0 );
	assert( row < getNbrOfRows() );

	assert( width >= 0 );
	assert( width < getNbrOfWidths() );


	int orderRow , orderColumn , orderWidth , orderHeight , tmp , m , n ;
	
	if( row > column ) {
		orderRow = column;
		orderColumn = row;
	} else {
		orderRow = row;
		orderColumn = column;
	}
	
	
	m =  orderColumn * (orderColumn + 1) / 2 + orderRow;

	if( width > height ) {
		orderWidth = height;
		orderHeight = width;
	} else {
		orderWidth = width;
		orderHeight = height;	
	}
		
	n =  orderHeight * (orderHeight + 1) / 2 + orderWidth;
		
	if( m > n ) {
		tmp = n;
		n = m;
		m = tmp;
	}
	
	return n * (n+1) / 2 + m;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method that reads a tensor in a file.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void tensorSymetric::readText(const string & fileName)
{
	ifstream fileStream;
	
	fileStream.open(fileName.c_str());
	
	if(!fileStream.is_open()) {
		cerr << "Error : in void tensorSymetric::readText(const string & fileName)." << endl;
		cerr << "Error : unable to open file with name " << fileName << endl;
		cerr << "Error : aborting reading now." << endl;
		clear();
		return;
	}
	
	readText(fileStream);
	
	if(fileStream.fail()) {
		cerr << "Error : in void tensorSymetric::readText(const string & fileName)." << endl;
		cerr << "Error : reading fails for file with name " << fileName << endl;
		clear();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that reads a tensor symetric in a stream.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void tensorSymetric::readText(istream & inStream) 
{
	
	int nbrOfRows;
	double value;
	
	inStream >> nbrOfRows;
	
	if(inStream.fail()) {
		cerr << "Error : in void tensorSymetric::readTxt(istream & inStream)." << endl;
		cerr << "Error : unable to read the number of rows in stream." << endl;
		cerr << "Error : aborting reading" << endl;
		clear();
		return;
	}
	
	setTensorSize(nbrOfRows);
	
	for(int i=0 ; i<getStorageSize() ; i++) {
		inStream >> value;
		
		if(inStream.fail()) {
			cerr << "Error : in void tensorSymetric::readTxt(istream & inStream)." << endl;
			cerr << "Error : unable to read data number " << i << "." << endl;
			cerr << "Error : aborting reading" << endl;
			clear();
			return;
		}
		
		setStorage(i,value);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that sets the size of the tensor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void tensorSymetric::setTensorSize(int nbrOfRows) 
{
	tensorNonSparse::setTensorSize(nbrOfRows , nbrOfRows , nbrOfRows , nbrOfRows);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that sets the size of the tensor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void tensorSymetric::setTensorSize(const int & nbrOfRows , const int & nbrOfColumns , const int & nbrOfWidths , const int & nbrOfHeights)
{
	assert( nbrOfRows == nbrOfColumns );
	assert( nbrOfWidths == nbrOfHeights );
	assert( nbrOfWidths == nbrOfRows );
	
	tensorNonSparse::setTensorSize(nbrOfRows , nbrOfColumns , nbrOfWidths , nbrOfHeights);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that writes down the tensor in a human readable way
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void tensorSymetric::writePretty(ostream & outStream , int width , int precision) const
{
	int i , j , k , l , nbrOfRows;
	
	nbrOfRows = getNbrOfRows();
	
	
	
	for(i=0 ; i < nbrOfRows ; i++) {
			
			for(j=0 ; j < nbrOfRows ; j++) {
				
				outStream << "Line " << setw(5) << i << " x" << setw(4) << j << " : " << setw(4) << i*nbrOfRows + j << " : ";
	
				if(j<i) {
					for(k=0 ; k < nbrOfRows ; k++) {		
						for(l=0 ; l < nbrOfRows ; l++) {		
							outStream << setw(width) <<"x";
						}
					} // fin du for k for l.
					
					outStream << endl;
					continue;
				}	// fin diu if i < j.

				for(k=0 ; k < nbrOfRows ; k++) {
				
					if( k < i) {
						for(l=0 ; l < nbrOfRows ; l++) {		
							outStream << setw(width) <<"x";
						}
						continue;
					} // fin du if k < i.
				
					
					for(l=0 ; l < nbrOfRows ; l++) {
					
						if( l < k ) {
							outStream << setw(width) <<"x";
							continue;
						}

						if(k==i && l<j) {
							outStream << setw(width) <<"x";
							continue;
						} 
																
						outStream << setw(width) <<  setprecision(precision) << getCoefficient(i,j,k,l);
						
					} // fin du for l.
				} // fin du for k.	
				
				
				outStream << endl;
			} // fin du for j.
		} // fin du for i.

/*	for(i=0 ; i < nbrOfRows ; i++) {
		for(j=0 ; j < nbrOfRows ; j++) {
			
			if(j<i) {
				for(k=0 ; k < nbrOfRows ; k++) {		
					for(l=0 ; l < nbrOfRows ; l++) {		
						outStream << setw(width) <<"x";
					}
				} // fin du for k for l.
				outStream << endl;
				continue;
			}	// fin diu if i < j.
			
			
			for(k=0 ; k < nbrOfRows ; k++) {
				
				if( k < i) {
					for(l=0 ; l < nbrOfRows ; l++) {		
						outStream << setw(width) <<"x";
					}
					continue;
				} // fin du if k < i.
				
				
				for(l=0 ; l < nbrOfRows ; l++) {
					
					if( l < k ) {
						outStream << setw(width) <<"x";
						continue;
					}
					
					if(k==i && l<j) {
						outStream << setw(width) <<"x";
						continue;
					} 
					
					outStream << setw(width) <<  setprecision(precision) << getCoefficient(i,j,k,l);
					
				} // fin du for l.
			} // fin du for k.	
			
			
			outStream << endl;
		} // fin du for j.
		} // fin du for i.
		*/
	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that writes the tensor with the only data that it contains.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void tensorSymetric::writeText(ostream & outStream , int precision) const
{
	outStream << getNbrOfRows() << endl << endl;
	
	for(int i =0 ; i < getStorageSize() ; i++) {
		outStream << setprecision(precision) << getStorage(i) << endl;
	}
}