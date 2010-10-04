/* 
* The martix library of the A.S.P.I.C. 
 * Written and directed by François Lodier francois.lodier@gmail.com.
 * Modified by Frederic Legoll nov06
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
#include "matrixFull.h"
#include <cstdlib>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The default constructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull::matrixFull(void)
: matrixNonSparse()
{
	setMatrixType(FULL);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The constructor with the sizes.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull::matrixFull(const int & nbrOfRows , const int & nbrOfColumns)
: matrixNonSparse()
{
	setMatrixType(FULL);
	setMatrixSize(nbrOfRows,nbrOfColumns);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The constructor with the sizes.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull::matrixFull(const matrixFull & mFull)
: matrixNonSparse()
{
	setMatrixType(FULL);
	copy(mFull);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The destructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull::~matrixFull(void)
{
	;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Method that computes the total storage size.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int matrixFull::computeStorageSize(void) const
{
	return getNbrOfRows() * getNbrOfColumns();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Method that computes the position in the storage array of the element row x column.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int matrixFull::hash(const int & row , const int & column) const
{
	assert( row >= 0 );
	assert( row < getNbrOfRows() );
	
	assert( column >= 0 );
	assert( column < getNbrOfColumns() );
	
	return row * getNbrOfColumns() + column;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that writes down the matrix in a human readable way
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void matrixFull::writePretty(ostream & outStream , int width , int precision) const
{
	int i , j , nbrOfColumns , nbrOfRows;	
	
	nbrOfRows = getNbrOfRows();
	nbrOfColumns = getNbrOfColumns();
	
	for(i=0 ; i < nbrOfRows ; i++) {
		
		for(j=0 ; j < nbrOfColumns ; j++) {
			
				outStream << setw(width) << setprecision(precision) << getCoefficient(i,j);
			
		}
		
		outStream << endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that computes the norm of a column of a full matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MatrixFullNormColumn(const matrixFull & A, int column)
{ 
double norm;
double norm2 = 0.0;
int i;
if (column>=A.getNbrOfColumns()){
   cerr << "Error in the norm of the column : the column does not exist." << endl;
   exit(1);
}
for (i=0;i<A.getNbrOfRows();i++){
norm2 += A(i,column) * A(i,column);   
}
norm = sqrt(norm2);
return norm;
} 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that computes the scalar product of two column of a full matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MatrixFullScalarProduct(const matrixFull & A, int i, int j)
{ double result = 0.0;
  int k;
  if (i>=A.getNbrOfColumns()){
   cerr << "Error in the scalar product of two columns : one column does not exist." << endl;
   exit(1);
}
if (j>=A.getNbrOfColumns()){
   cerr << "Error in the scalar product of two columns : one column does not exist." << endl;
   exit(1);
}
  for (k=0;k<A.getNbrOfRows();k++){
     result += A(k,i) * A(k,j);
  }
  return result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that takes the transpose of a full matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull matrixFull::TransposeMatrix(void) const
{ int i, j;
  matrixFull TA(getNbrOfColumns(),getNbrOfRows());
  
  for (i=0;i<getNbrOfColumns();i++){
     for (j=0;j<getNbrOfRows();j++){
	    TA(i,j) = getCoefficient(j,i);
	 }
  }  
  
  return TA;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that takes the signature of an element of a full matrix 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int matrixFull::MatrixSignElement(int i, int j) const
{ int k, test, sign;

  if ((i<0)||(i>=getNbrOfRows())||(j<0)||(j>=getNbrOfColumns())){
  cerr << "Error : one of the arguments of MatrixSignElement is out of the dimension of the matrix." << endl;
  exit (1);
  }
  k = i + j;
  test = k%2;
  if (test == 1){
     sign = -1;
  } 
  else {
  sign = 1;
  }

  return sign;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that find the minor of a square full matrix 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void matrixFull::MinorForSquareMatrix(matrixFull & M1, int line, int column) const
{ int i, j, a, b;

   assert (getNbrOfRows() == getNbrOfColumns());
   if ((line<0)||(line>=getNbrOfRows())||(column<0)||(column>=getNbrOfRows())){
      cerr << "Error : one of the arguments of MinorForSquareMatrix is out of the dimension of the matrix." << endl;
      exit (1);
  }
  
  matrixFull M(getNbrOfRows(),getNbrOfRows());
  //matrixFull M1(getNbrOfRows()-1,getNbrOfRows()-1);
  
  a = 0;
  for (i=0;i<getNbrOfRows();i++){
     b = 0;
	 if (i!=line){
	    for (j=0;j<getNbrOfRows();j++){
		   if (j!=column){
		      M(a,b) = getCoefficient(i,j);
			  b++;
		   }
		}
		a++;
	 }
  }

  for (i=0;i<getNbrOfRows()-1;i++){
     for (j=0;j<getNbrOfRows()-1;j++){
	    M1(i,j) = M(i,j);
	 }
  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that computes the determinant of a full matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double matrixFull::Determinant(void) const
{  
   assert (getNbrOfRows() == getNbrOfColumns());
   
   int k;
   double det = 0.0;
   if (getNbrOfRows()==1){ 
      det = getCoefficient(0,0);
   }	  
   else {
      int dim = getNbrOfRows() - 1;
	  matrixFull Aux(dim,dim);
      for (k=0;k<getNbrOfRows();k++){
	     //for (i=0;i<dim;i++){
		    //for (j=0;j<dim;j++){
			   MinorForSquareMatrix(Aux,0,k);
			//}// end of j loop
		 //}// end of i loop
		 det += MatrixSignElement(0,k) * getCoefficient(0,k) * Aux.Determinant();
	  }// end of k loop
   }// end else
  
   return det;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that computes the comatrix of a full matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull matrixFull::Comatrix(void) const
{  
   assert (getNbrOfRows() == getNbrOfColumns());
   int i, j;
   matrixFull C(getNbrOfRows(),getNbrOfRows());
   matrixFull M(getNbrOfRows(),getNbrOfRows());
   
   
   for (i=0;i<getNbrOfRows();i++){
      for (j=0;j<getNbrOfRows();j++){
	     MinorForSquareMatrix(M,i,j);
	     C(i,j) = MatrixSignElement(i,j) * M.Determinant();
	  }
   }
   return C;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that computes the inverse of a full matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull matrixFull::Inverse(void) const
{ assert (getNbrOfRows() == getNbrOfColumns());

  matrixFull inv(getNbrOfRows(),getNbrOfRows());
  matrixFull B(getNbrOfRows(),getNbrOfRows());
  int i,j;
  
  double det;
  det = Determinant();
  if (det==0){
     cerr << "The matrix is not inversible." << endl;
	 exit (1);
  }
  else{
     for (i=0;i<getNbrOfRows();i++){
	    for (j=0;j<getNbrOfRows();j++){
		   inv(i,j) = Comatrix().TransposeMatrix()(i,j) / det;
		}// end of j loop
	 }// end of i loop
  }// end else
  
  return inv;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The method that computes the LU decomposition of a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void matrixFull::LUdecomposition(matrixFull & LU, int *indx, double d) const
{ 
  assert (getNbrOfRows() == getNbrOfColumns());
  assert (LU.getNbrOfRows() == getNbrOfRows());
  assert (LU.getNbrOfColumns() == getNbrOfRows());
  
  const double TINY = 1.0e-20;
  int i, imax, j, k;
  double big, dum, sum, temp;
  double *vv;
  
  cout <<"Performing the LU decomposition of the matrix." << endl;
  
  
  vv = new double [getNbrOfRows()];
  //indx = new int [getNbrOfRows()];
  
  for (i=0;i<getNbrOfRows();i++){
     for (j=0;j<getNbrOfRows();j++){
	    LU(i,j) = getCoefficient(i,j);
	 }
  }
  d = 1.0;  
  for (i=0;i<getNbrOfRows();i++){
     big = 0.0;
	 for (j=0;j<getNbrOfRows();j++)
	   // ALH - more brackets around "=" to avoid g++ 4.1 warnings.
	    if ((temp=fabs(LU(i,j)))>big) big = temp;
		if (big==0.0){
		   cerr << "Error : Singular matrix in the method LUdecomposition." << endl;
		   exit (1);
		}// end if
	 //}//end of j loop 
	 vv[i] = 1.0 / big;
  }// end of i loop
  for (j=0;j<getNbrOfRows();j++){
     for (i=0;i<j;i++){
	    sum = LU(i,j);
		for (k=0;k<i;k++){
		   sum -= LU(i,k) * LU(k,j);
		}// end of k loop
		LU(i,j) = sum;
	 }// end of i loop
	 big = 0.0;
	 for (i=j;i<getNbrOfRows();i++){
        sum = LU(i,j);
		for (k=0;k<j;k++){
		   sum -= LU(i,k) * LU(k,j);
		}// end of k loop
		LU(i,j) = sum;
		if ((dum=vv[i]*fabs(sum))>=big){
		   big = dum;
		   imax = i;
		}// end if   
	 }// end of i loop
	 if (j != imax){
	    for (k=0;k<getNbrOfRows();k++){
		   dum = LU(imax,k);
		   LU(imax,k) = LU(j,k);
		   LU(j,k) = dum;
		}// end of k loop
		d = - d;
		vv[imax] = vv[j];
	 }// end if
	 indx[j] = imax;
	 if (LU(j,j)==0.0) LU(j,j) = TINY;
	 if (j != getNbrOfRows()-1){
	    dum = 1.0 / (LU(j,j));
		for (i=j+1;i<getNbrOfRows();i++){
		   LU(i,j) *= dum;
		}// end of i loop
	 }// end if
  }// end of j loop

delete [] vv;

cout <<"The LU decomposition is done" << endl;

}

// Function added by F. Legoll
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Affectation operator
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrixFull & matrixFull::operator= (const matrixFull & mFull)
{
	copy(mFull);
	return *this;
}
