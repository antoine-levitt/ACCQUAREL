/* 
 * The chemics library of A.S.P.I.C. 
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
#ifndef _GAUSSIAN_BASIS_SHELL_
#define _GAUSSIAN_BASIS_SHELL_

#include "contractions.h"
#include <dpoint.h>
#include "gaussianBasisFunction.h"
#include "shellType.h"
#include <iostream>
using namespace std;


/**
 * Class that manage shell.
 */
class gaussianBasisShell : public contractions , public shellType
{
private:
	
	/**
	 * Le centre de la couche.
	 */
	dpoint<3> Center;

protected:

	/**
	 * Method that performs copy.
	 */
	void copy(const gaussianBasisShell & gbs);
		
public:

	/**
   * Default constructor.
	 */
	gaussianBasisShell(void);

	/**
	 * Copy constructor.
	 */
	gaussianBasisShell(const gaussianBasisShell & gbs);

	/**
   * Destructor.
	 */
	virtual ~gaussianBasisShell(void);
	
	
	/**
	 * The clear method.
	 */
	void clear(void);
	
	/**
	 * Method GET for a basis function.
	 *
	 * As a shell is a set of basis functions, this method allow
	 * user to extract on of them.
	 *
	 * @warning the parameter item must be in the good range.
	 */
	gaussianBasisFunction getBasisFunction(int item) const;

	/**
	 * Method GET for the center of the shell.
	 */
	const dpoint<3> & getCenter(void) const;

	/**
	 * Method SET for the center of the shell.
	 */
	void setCenter(const dpoint<3> & center);
	
	/**
	 * Method WRITE.
	 */
	void write(ostream & out) const;
};

/**
 * Extern operator << .
 */
ostream & operator << (ostream & out , const gaussianBasisShell & gbs);

#endif

