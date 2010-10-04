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
#ifndef _BASIS_ELEMENT_
#define _BASIS_ELEMENT_

#include <containor.h>
#include "contractions.h"
#include "shellType.h"

/**
 * Class to represent an atom in a gaussian basis.
 */
class basisElement 
{
	
public :
	
	/**
	 * Enumération pour le type de base manipulé.
	 * En fait, aujourd'hui tout les éléments sont de type Gaussian 
	 * mais qui nous dit qu'on ne pourra jamais faire de Slater opu autre ...
	 */
	typedef enum _basisType {
		GAUSSIAN , 
		SLATER
	} basisType;
	
private:
	
 /**
	* This  string contains the name of the  gaussian basis in 
	* wich the element is described.
	*/
	string BasisName;

	/**
	 * This string contains the key of the chemical element
	 * that is described in the gaussian basis.
	 */
	string BasisElementKey;

	/**
	 * This contains information about the shells
	 * that contains the gaussian basis element.
	 */
	containor<contractions> Contractions4Shells;
	
	/**
   * This contains iformations about the shells Types
	 * that are contained in the gaussian basis element.
	 */
	containor<shellType> ShellTypes;

	/**
	 * Ici on stocke le type de base que l'on est entraint de manipuler.
	 */
	basisType Type;
	
protected:
	
	/**
	 * Method that copy a gaussianBasisElement.
	 */
	void copy(const basisElement & gBasisElement);

public:

	/**
	 * Le constructeur.
	 */
	basisElement(void);

	/**
	 * Le destructeur.
	 */
	~basisElement(void);

	/**
   * Method that clears the object.
	 */
	void clear(void);
	
	/**
	 * Method to know of an element is empty.
	 */
	bool empty(void) const;

	/**
	 * Accesseur au nom de la base.
	 */
	const string & getBasisName(void) const;

	/**
	 * Method GET for the chemical element key.
	 */
	const string & getBasisElementKey(void) const;

	/**
	 * Method GET for a shell contractions.
	 */
	const contractions & getContractions4Shell(int item) const;
	
	/**
	 * Method GET for all the shells containors.
	 */
	const containor<contractions> & getContractions4Shells(void) const;
	
	/**
	 * Method GET for the number of contractions of a given Shell.
	 */
	int getNbrOfContractions4Shell(int item) const;

	/**
	 * Method GET for the number of basis functions of a shell 
	 * in the element.
	 */
	int getNbrOfBasisFunctions4Shell(int item) const;

	/**
	 * Method GET for the number of shells.
	 */
	int getNbrOfShells(void) const;

	/**
   * Méthode pour accéder au type de couche de la couche item.
	 */
	const shellType & getShellType(int item) const;

	/**
   * Méthode pour connaitre le nom de la clé (de type de couche) de la couche item.
	 */
	const string & getShellTypeKey(int item) const;

	/**
	 * Method GET for the polynomes degree of a given shell and basis function.
	 */
	const ipoint<3> & getMonomeDegree(int shell , int function) const;

	/**
	 * Method GET for all the shell types.
	 */
	const containor<shellType> & getShellTypes(void) const;

	/**
	 * Méthode pour savoir à quelle type de base nous avons à faire.
	 */
	const basisType & getBasisType(void) const;
	
	/**
   * Operator=.
	 */
	basisElement & operator=(const basisElement & element); 

	/**
   * Method set for the name of the basis.
	 */
	void setBasisName(const string & basisName);

	/**
	 * Method SET for the chemical element key.
	 */
	void setBasisElementKey(const string & basisElementKey);
		
	/**
   * Method SET for a basis shell contraction.
	 */
	void setContractions4Shell (int item , const contractions & shellContractions);
	
	/**
   * Method SET for all shell contractions.
	 * 
	 * @warning as the number of contractions must be the same as the 
	 * number of shell type keys this method does not performs
	 * memory allocation. Call setNbrOfShells() first ! 
	 */
	void setContractions4Shells (const containor<contractions> & shellContractions);

	/**
	 * Method SET for the whole element.
	 *
	 * This method construct the object perfonrming a connection to the data base.
	 */
	void setGaussianBasisElement(const string & basisName , const string & basisElementKey);
	
	/**
	 * Method SET for the whole element.
	 *
	 * This method construct the object performing copy of the argument.
	 */
	void setBasisElement(const basisElement & element);	
	
	/**
	 * Method SET for the number of shells.
	 */
	void setNbrOfShells(int nbrOfShells);
	
	/**
	 * Method SET for a shell type.
	 */
	void setShellType(int item , const string & shellTypeKey);

	/**
	 * Method SET for all the shell type keys.
	 * 
	 * @warning as the number of contractions must be the same as the 
	 * number of shell type keys this method does not performs
	 * memory allocation. Call setNbrOfShells() first ! 
	 */
	void setShellTypes(const containor<shellType> & shellTypeKeys);

	/**
   * Method to write thet gaussianBasisElement.
	 */
	void write(ostream & out) const;
};

/**
 * operator << for a basisElement
 */
ostream & operator<<(ostream & out , const basisElement & element);
#endif

