/* 
 * The containor library of the A.S.P.I.C. 
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
#ifndef _CONTAINOR_2D_
#define _CONTAINOR_2D_

#include "containor.h"
#include "ipoint.h"
#include <iostream>
#include <iomanip>
using namespace std;

/**
 * Classe pour manipuler des tableaux bi dimensionnels.
 */
template<typename Type>
class containor2D : public containor<Type> {
private:
	
	/**
	 * Le nombre de lignes du tableau.
	 */
	int NbrOfRows; 	
	
	/**
	 * Le nombre de colones du tableau.
	 */
	int NbrOfColumns; protected: 	
	
	/**
	 * Méthode qui permet de faire la copie.
	 *
	 * @param c le tableau à copier dans l'objet appelant.
	 */
	inline void copy(const containor2D<Type> & c)
	{
		NbrOfRows = c.getNbrOfRows();
		NbrOfColumns = c.getNbrOfColumns();
		containor<Type>::copy(c);
	} 	
	
	/**
	 * Méthode qui permet de faire le hash.
	 *
	 * @param row la ligne de l'élément du tableau auquel on souhaite accéder.
	 * @param column la colone de l'élément du tableau auquel on souhaite accéder.
	 *
	 * @return la position le l'élément (row , column) dans le containor qui stoque 
	 * les valeurs des éléments du tableau.
	 */ 	
	inline int hash(const int & row , const int & column) const
	{
		assert(row >= 0);
		assert(row < getNbrOfRows());
		assert(column >= 0);
		assert(column < getNbrOfColumns()); 		
		return row + getNbrOfRows() * column;
	} 

public: 	
	
	/**
	 * Constructeur par défaut.
	 *
	 * Construit un objet totalement vide qui ne peut pas 
	 * contenir de données.
	 */
	inline containor2D(void)
		: containor<Type>() ,
		NbrOfColumns(0) ,
		NbrOfRows(0)
	{
		;
	} 	
	
	/**
	 * Constructeur avec spécification des tailles.
	 *
	 * @param nbrOfRows le nombre de lignes que le tableau pourra contenir.
	 * @param nbrOfColumns le nombre de colones que le tableau pourra contenir.
	 */
	inline containor2D(const int & nbrOfRows , const int & nbrOfColumns)
		: containor<Type>() ,
		NbrOfColumns(0) ,
		NbrOfRows(0)
	{
		this->setSize(nbrOfRows,nbrOfColumns);
	} 	
	
	/**
	 * Constructeur de copie.
	 *
	 * @param c le tableau à copier.
	 */
	inline containor2D(const containor2D<Type> & c)
		: containor<Type>() ,
		NbrOfColumns(0) ,
		NbrOfRows(0)
	{
		copy(c);
	} 	
	
	/**
	 * Destructeur.
	 */
	inline virtual ~containor2D(void)
	{
		clear();
	} 	
	
	/**
	 * Méthode qui nettoie totalement l'objet.
	 * 
	 * Cette méthode permet de reinitialiser totalement le tableau 
	 * à zéro. Après l'appel de celle ci, plus aucune donnée ne peut
	 * etre stoquée dans le tableau.
	 */ 	inline void clear(void)
	{
		NbrOfColumns = 0;
		NbrOfRows = 0;
		containor<Type>::clear();
	} 	
	
	/**
	 * Méthode pour accéder aux élement du tableau.
	 *
	 * @param row la ligne de l'élément dont on souhaite connaitre la valeur.
	 * @param column la colone de l'élément dont on souhaite connaitre la valeur.
	 *
	 * @warning il faut que row soit non négatif et strictement inférieur aux
	 * nombre de lignes que peut contenir l'objet.
	 *
	 * @warning il faut que column soit non négatif et strictement inférieur aux
	 * nombre de colones que peut contenir l'objet.
	 */
	inline Type & getData(const int & row , const int & column) 
	{
		assert(row >=0);
		assert(column >=0); 		assert(row < getNbrOfRows());
		assert(column < getNbrOfColumns()); 		return containor<Type>::getData(hash(row,column));
	} 	
	
	/**
	 * Méthode constante pour accéder aux éléments du tableau.
	 *
	 * @param row la ligne de l'élément dont on souhaite connaitre la valeur.
	 * @param column la colone de l'élément dont on souhaite connaitre la valeur.
	 *
	 * @warning il faut que row soit non négatif et strictement inférieur aux
	 * nombre de lignes que peut contenir l'objet.
	 *
	 * @warning il faut que column soit non négatif et strictement inférieur aux
	 * nombre de colones que peut contenir l'objet.
	 */
	inline const Type & getData(const int & row , const int & column) const
	{
		assert(row >=0);
		assert(column >=0); 		assert(row < getNbrOfRows());
		assert(column < getNbrOfColumns()); 		return containor<Type>::getData(hash(row,column));
	} 	
	
	/**
	 * Méthode pour connaitre le nombre de lignes du tableau.
	 *
	 * @return le nombre de lignes dans le tableau.
	 */
	inline const int & getNbrOfRows(void) const
	{
		return NbrOfRows;
	} 	
	
	/**
	 * Méthode pour connaitre le nombre de colones du tableau.
	 *
	 * @return le nombre de colones dans le tableau.
	 */
	inline const int & getNbrOfColumns(void) const
	{
		return NbrOfColumns;
	} 	
	
	/**
	 * Méthode pour connaitre la taille du tableau.
	 *
	 * @return une paire d'entiers qui contient le nombre de lignes
	 * respectivement le nombre de colones) du tableau pour la première composante 
	 * (respectivement pour la seconde).
	 */
	inline ipoint<2> & getSizes(void) const
	{
		ipoint<2> sizes(getNbrOfRows() , getNbrOfColumns());
		return sizes;
	} 	
	
	/**
	 * Méthode pour connaitre la taille du tableau dans une direction donnée.
	 *
	 * @param dimension la dimension dans laquelle on souhaite connaitre la taille.
	 *
	 * @return le nombre de lignes ou de colones dans le tableau suivant le numéro de
	 * la dimension demandée.
	 *
	 * @warning la dimension doit prendre la valeur 0 ( pour accéder aux nombre de lignes)
	 * ou 1 (pour accéder au nombre de colones) mais ne peut prendre aucune autre valeur.
	 */
	inline const int & getSize(int dimension) const
	{
		assert(dimension >= 0);
		assert(dimension < 2);

		if(dimension == 0) {
			return getNbrOfRows();
		} else {
			return getNbrOfColumns();
		}		
	} 	
	
	/**
	 * Opérateur d'affectation.
	 *
	 * @param c le containor à copier.
	 */
	inline containor2D<Type> & operator=(const containor2D<Type> & c)
	{
		copy(c);
		return *this;
	} 	
	
	/**
	 * Opérateur constant pour accéder aux éléments du tableau.
	 *
	 * @param row la ligne de l'élément dont on souhaite connaitre la valeur.
	 * @param column la colone de l'élément dont on souhaite connaitre la valeur.
	 *
	 * @warning il faut que row soit non négatif et strictement inférieur aux
	 * nombre de lignes que peut contenir l'objet.
	 *
	 * @warning il faut que column soit non négatif et strictement inférieur aux
	 * nombre de colones que peut contenir l'objet.
	 */
	inline const Type & operator()(const int & row ,const int & column) const
	{
		assert(row >=0);
		assert(column >=0); 		assert(row < getNbrOfRows());
		assert(column < getNbrOfColumns()); 		return getData(row,column);
	} 	
	
	/**
	 * Opérateur d'accès aux éléments du tableau.
	 *
	 * @param row la ligne de l'élément dont on souhaite connaitre la valeur.
	 * @param column la colone de l'élément dont on souhaite connaitre la valeur.
	 *
	 * @warning il faut que row soit non négatif et strictement inférieur aux
	 * nombre de lignes que peut contenir l'objet.
	 *
	 * @warning il faut que column soit non négatif et strictement inférieur aux
	 * nombre de colones que peut contenir l'objet.
	 */
	inline Type & operator()(const int & row , const int & column)
	{
		assert(row >=0);
		assert(column >=0); 		assert(row < getNbrOfRows());
		assert(column < getNbrOfColumns()); 		return getData(row,column);
	} 

	/**
	 * Méthode pour fixer les valeurs du tableau à partir d'un autre tableau.
	 *
	 * @param array le tableau que l'on souhaite recopier dans l'objet.
	 */
	inline void setDatas(const containor2D<Type> & array)
	{
		copy(array);
	}

	/**
	 * Méthode pour remplir le conteneur avec une seule valeure.
	 *
	 * @param value la valeur avec laquelle on souhaite remplir le containor.
	 */
	inline void setDatas(const Type & value)
	{
		containor<Type>::setDatas(value);
	}
	
	/**
	 * Méthode pour spécifier la taille du tableau
	 *
	 * @param nbrOfRows le nombre de lignes pour le tableau.
	 * @param nbrOfColumns le nombre de colones pour le tableau.
	 *
	 * @warning il faut que le paramêtre nbrOfLines soit strictement
	 * positif.
	 *
	 * @warning il faut que le parametre nbrOfColumns soit strictement 
	 * positif.
	 *
	 * @see clear(void)
	 */
	inline void setSizes(const int & nbrOfRows , const int & nbrOfColumns)
	{
		assert(nbrOfRows > 0);
		assert(nbrOfColumns > 0); 		
		NbrOfRows = nbrOfRows;
		NbrOfColumns = nbrOfColumns; 		
		containor<Type>::setSizes(nbrOfRows * nbrOfColumns);
	} 	
	
	/**
	 * Méthode qui permet d'afficher le conteneur sous une forme
	 * lisible pour l'homme.
	 *  
	 * @param out le flux dans lequel on veut écrire le contrainor.
	 */
	inline virtual void writePretty(ostream & out) const
	{
		int row , column , nbrOfRows , nbrOfColumns;
		
		int width = 8;

		nbrOfRows = getNbrOfRows();
		nbrOfColumns = getNbrOfColumns();

		for( row = 0 ; row < nbrOfRows ; row++) {
			
			for(column = 0 ; column < nbrOfColumns ; column++) {	
				out << setw(width) << getData(row,column);
			}

			if(row != (nbrOfRows - 1)) 
				out << endl;
		}
	}

};

#endif  
