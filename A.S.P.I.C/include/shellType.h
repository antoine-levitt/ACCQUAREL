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
#ifndef _GAUSSIAN_BASIS_SHELL_TYPE_
#define _GAUSSIAN_BASIS_SHELL_TYPE_


#include <containor.h>
#include <ipoint.h>

class shellType
{
private:
	
	/**
	 * La chaine de charactère pour qui contient la clé du type de couche dans
	 * la base de données.
	 */
	string ShellTypeKey;
	
	/**
   * La chaine de charactère qui contient le nom de la couche.
	 */
	string ShellTypeName;

	/**
   * Les degres des monomes qui composent la couche.
	 */
	containor< ipoint<3> > MonomeDegrees;
	
protected:

	/**
	 * Méthode pour faire la copie d'un type de couche.
	 */
	void copy(const shellType & shell);

public:

	/**
   * Constructeur par défaut.
	 */
	shellType(void);

	/**
	 * Constructeur de copie.
	 *
	 * @param shell le type de couche avec lequel l'objet nouvellement créer sera instancié.
	 */
	shellType(const shellType & shell);

	/**
	 * Destructeur.
	 */
	~shellType(void);

	/**
	 * Méthode pour remmettre le type de couche à zéro.
	 */
	void clear(void);
	
	/**
	 * Méthode pour savoir si le type de couche est vide.
	 */
	bool empty(void) const;
	
	/**
	 * Méthode pour connaitre le plus haut degré de la couche.
	 */
	int getHigherDegree(void) const;
	
	/**
	 * Méthode pour connaitre le degré du monome item de la couche.
	 *
	 * @warning il faut que le parametre item soit compris entre
	 * 0 et le nombre de monomes présents dans la couche.
	 */
	const ipoint<3> & getMonomeDegree(int item) const;
	
	/**
	 * Méthode pour accéder aux degrés de la couche comme un tableau
	 * de degrés.
	 */
	const containor< ipoint<3> > & getMonomeDegrees(void) const;
	
	/**
	 * Méthode pour connaitre le nombre de fonctions de base dans la couche.
	 * Cette méthode revient à compter le nombre de degré de monomes présents dans la 
	 * couche.
	 */
	int getNbrOfBasisFunctions(void) const;
	
	/**
	 * Méthode pour connaitre le nom de la clé du type 
	 * de couche dans la base de données. Si l'objet n'a pas
	 * été crée par la lecture dans la base de données alors
	 * la cahine retournée est vide.
	 */
	const string & getShellTypeKey(void) const;
	
	/**
	 * Méthode pour connaitre le nom du type de manipulé couche.
	 */
	const string & getShellTypeName(void) const;

	/**
	 * Opérateur d'affectation.
	 */
	shellType & operator=(const shellType & shell);
	
	/**
	 * Méthode pour spécifier les degrés de tout les monomes.
	 *
	 * Cette méthode recopie le tableau passés en argument comme étant
	 * l'ensemble de degrés des monomes qui composent la couche.
	 */
	void setMonomeDegrees(const containor< ipoint<3> > & monomeDegrees);

	/**
	 * Méthode pour modifier le degré d'un monome.
	 *
	 * @param item le numéro du monome dont on souhaite modifier le degré.
	 *
	 * @param monomeDegree la valeur du degré à affecter au monome item. 
	 */
	void setMonomeDegree(int item , const ipoint<3> & monomeDegree);

	/**
	 * Méthode pour modifier le degré d'un monome.
	 *
	 * @param item le numéro du monome dont on souhaite modifier le degré.
	 *
	 * @param monomeDegreeX la valeur du degré à affecter au monome item pour la composante en X. 
	 * @param monomeDegreeY la valeur du degré à affecter au monome item pour la composante en Y. 
	 * @param monomeDegreeZ la valeur du degré à affecter au monome item pour la composante en Z. 
	 */
	void setMonomeDegree(int item , int monomeDegreeX , int monomeDegreeY , int monomeDegreeZ);

	/**
   * Méthode pour modifier le nombre de fonctions de bases présentent dans la couche.
	 * 
	 * @param nbrOfBasisFunctions le nombre de fonctions de base présentes dans la couche.
	 */
	void setNbrOfBasisFunctions(int nbrOfBasisFunctions);

	/**
	 * Méthode pour retrouver un type de couche à partir de sa clé dans la base de données.
	 *
	 * @parma shellTypeKey la clé à rechercher dans la base de données.
	 */
	void setShellType(const string & shellTypeKey);

	/**
	 * Méthode pour recopier les valeurs d'un type de couche dans l'objet.
	 *
	 * @param shell le type de couche à recopier dans l'objet.
	 */
	void setShellType(const shellType & shell);

	/**
	 * Méthode pour modifier la clé du type de couche.
	 *
	 * @param shellTypeKey la clé à attribuer au type de couche.
	 */
	void setShellTypeKey(const string & shellTypeKey);
	
	/**
	 * Méthode pour modifier le nom du type de couche.
	 *
	 * @param shellTypeName le nom a attribuer au type de couche.
	 */
	void setShellTypeName(const string & shellTypeName);
	
	/**
	 * Méthode pour écrire le type de couche dans un format humain.
	 */
	void writeHuman(ostream & out , const string & linePrefix = "") const;
};


/**
 * Opérateur externe pour l'écriture dans un flux.
 *
 * @param outStream le flux dans lequel on souahite écrire la molécule.
 *
 * @param shell le type de couche à écrire dans le flux.
 */
ostream & operator<<(ostream & outStream , const shellType & shell);

#endif

