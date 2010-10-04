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
#ifndef _containorTri2D_
#define _containorTri2D_

#include "containor.h"
#include "ipoint.h"
#include <iomanip>
#include <iostream>

using namespace std;

/**
 * Classe permettant de contenir des tableaux bi dimensionnels
 * triangulaires.
 *
 * En fait cela revient a dire que nous souhaitons stoquer des choses 
 * de la forme suivante :
 *
 * Element0-0     x             x             x           ...
 * Element1-0     Element1-1    x             x           ...
 * Element2-0     Element2-1    Element2-2    x           ...
 * Element3-0     Element3-1    Element3-2    Element3-3  x      ...
 * ...
 *
 * Par convention nous choisissons de stocker la partie inférieure
 * (ie row >= column), même si nous aurions pu faire l'inverse.
 *
 * @warning Comme nous n'avons pas de notion de symétrie à ce stade, les méthodes
 * de la classe comportent des contraintes fortes sur le fait que nous 
 * stoquions uniquement la partie inférieure : cela veut dire que les
 * tentative d'accés à la partie non définie du tableau (ie les éléments pour lesquels 
 * row < column) provoquerons des erreurs (et l'interruption du programe en mode de 
 * debuguage).
 */
template<typename Type>
class containorTri2D : public containor<Type>
{
private:	

	/**
	 * Variable contenant le nombre de lignes du tableau.
	 * Comme nous travaillons avec un conteneur triangulaire,
	 * le nombre de colones est forcément identique aux nombres
	 * de lignes.
	 */
	int NbrOfRows;

protected:	

	/**
	 * Méthode de copie.
	 * 
	 * @param c le containor à copier dans l'objet courrant.
	 */
	inline void copy(const containorTri2D<Type> & c)
	{
		NbrOfRows = c.getNbrOfRows();
		containor<Type>::copy(c);
	}

/**
	 * Méthode qui permet de retrouver la position d'un élément du tableau
	 * bi dimensionel dans le vecteur de stockage.
	 *
	 * @param row la ligne de l'élément.
	 *
	 * @param column la colone de l'élément.
	 *
	 * @return la position de cet élément dans le tableau qui sert à stoquer les
	 * valeurs des éléments.
	 *
	 * @warning il faut que le paramètre row soit non négatif et strictement
	 * inférieur au nombre de lignes de l'objet.
	 *
	 * @warning il faut que le paramètre column soit non négatif et inférieur ou 
	 * égal à la valeur de row puisque nous sommes dans un cas triagulaire.
	 *
	 */	
	inline int hash(const int & row , const int & column) const 
	{
		assert(row >= 0);
		assert(row < getNbrOfRows());

		assert(column >= 0);
		assert(column <=  row);

		return (row * (row+1)) / 2 + column;
	}

public:

	/**
	 * Constructeur par défaut.
	 *
	 * Construit un tableau totalement vide.
	 */
	inline containorTri2D(void)
		: containor<Type>() ,
		NbrOfRows(0)
	{
		;
	}

	/**
	 * Constructeur avec la taille.
	 * 
	 * @param nbrOfRows le nombre de lignes (et par conséquent de colones)
	 * que pourra contenir le tableau bi-dimensionnel.
	 *
	 * @warning le nombre de lignes doit etre strictement positif.
	 */
	inline containorTri2D(const int & nbrOfRows)
		: containor<Type>() ,
		NbrOfRows(0)
	{
		assert(nbrOfRows > 0);
		setSizes(nbrOfRows);
	}

	/**
	 * Constructeur de copie.
	 * 
	 * @param array le tableau qui devra initialiser l'objet crée.
	 */
	inline containorTri2D(const containorTri2D<Type> & array)
		: containor<Type>() ,
		NbrOfRows(0)
	{
		copy(array);
	}
	
	/**
	 * Destructeur.
	 */
	inline virtual ~containorTri2D(void)
	{
		clear();
	}

	/**
	 * Méthode qui permet de nettoyer le tableau.
	 * 
	 * @warning Après l'appel de cette méthode plus rien ne peut
	 * etre stocker dans le tableau, pour allouer un nouvel espace
	 * il faut utiliser la méthode void set_nbr_of_lines(int)
	 */
	inline void clear(void)
	{
		NbrOfRows = 0;
		containor<Type>::clear();
	}
	
	/**
	 * Méthode pour accéder aux éléments du tableau.
	 *
	 * @param row la ligne de l'élément recherché.
	 *
	 * @param column la colone de l'élément recherché.
	 *
	 * @return une référence vers l'élément (row , column).
	 *
	 * @warning la ligne de l'élément recherché doit etre non négative
	 * et strictement inférieure aux nombre de lignes que contient le
	 * tableau.
	 *
	 * @warning la colone de l'élément recherché doit etre non négative
	 * et inférieure ou égale à la ligne de celui ci.
	 */
	inline Type & getData(const int & row ,const int & column) 
	{
		assert(row >=0);
		assert(column >=0);

		assert(row < getNbrOfRows());
		assert(column <= row);

		return containor<Type>::getData(hash(row,column));
	}

	/**
	 * Méthode constante pour accéder aux éléments du tableau.
	 *
	 * @param row la ligne de l'élément recherché.
	 *
	 * @param column la colone de l'élément recherché.
	 *
	 * @return une référence constante vers l'élément (row , column).
	 *
	 * @warning la ligne de l'élément recherché doit etre non négative
	 * et strictement inférieure aux nombre de lignes que contient le
	 * tableau.
	 *
	 * @warning la colone de l'élément recherché doit etre non négative
	 * et inférieure ou égale à la ligne de celui ci.
	 */
	inline const Type & getData(const int & row ,const int & column) const
	{
		assert(row >=0);
		assert(column >=0);

		assert(row < getNbrOfRows());
		assert(column <= row);

		return containor<Type>::getData(hash(row,column));
	}

	/**
	 * Méthode pour connaitre le nombre de lignes et par conséquent 
	 * de colones que peut contenir le tableau.
	 * 
	 * @return le nombre de lignes (aussi le nombre de colones) pour lequel l'espace mémoire à
	 * été alloué. 
	 */
	inline const int & getNbrOfRows(void) const
	{
		return NbrOfRows;
	}

	/**
	 * Méthode pour connaitre le nombre de colones et par conséquent 
	 * de lignes que peut contenir le tableau.
	 * 
	 * @return le nombre de colones (aussi le nombre de lignes) pour lequel l'espace mémoire à
	 * été alloué. 
	 */
	inline const int & getNbrOfColumns(void) const
	{
		return NbrOfRows;
	}

	/**
	 * Méthode pour connaitre la taille du tableau
	 * 
	 * @return une paire d'entiers qui contient 
	 * le nombre de lignes et le nombre de colones  que contient le tableau.
	 */
	inline ipoint<2> getSizes(void) const
	{
		ipoint<2> sizes(getNbrOfRows());
		return sizes;
	}

	/**
	 * Méthode pour connaitre la taille du tableau
	 * 
	 * @return une paire d'entiers qui contient 
	 * le nombre de lignes et le nombre de colones  que contient le tableau.
	 */
	inline ipoint<2> getSize(const int & dimension) const
	{
		assert( dimension >= 0);
		assert( dimension < 2);
	
		return getNbrOfRows();
	}

	/**
	 * Opératuer pour l'affectation.
	 * 
	 * @param array le tableau à copier.
	 *
	 * @return une référence vers l'objet modifié.
	 */
	inline containorTri2D<Type> & operator= (const containorTri2D<Type> & c)
	{
		copy(c);
		return *this;
	}

	/**
	 * Opérateur constant pour accéder aux éléments du tableau.
	 *
	 * @param row la ligne de l'élément recherché.
	 *
	 * @param column la colone de l'élément recherché.
	 *
	 * @return une référence constante vers l'élément (row , column).
	 *
	 * @warning la ligne de l'élément recherché doit etre non négative
	 * et strictement inférieure aux nombre de lignes que contient le
	 * tableau.
	 *
	 * @warning la colone de l'élément recherché doit etre non négative
	 * et inférieure ou égale à la ligne de celui ci.
	 */
	inline const Type & operator()(const int &  row , const int & column) const
	{
		assert(row >= 0);
		assert(column >= 0);
		
		assert(row < getNbrOfRows());
		assert(column <= row);

		return getData(row,column);
	}

	/**
	 * Opérateur pour accéder aux éléments du tableau.
	 *
	 * @param row la ligne de l'élément recherché.
	 *
	 * @param column la colone de l'élément recherché.
	 *
	 * @return une référence vers l'élément (row , column).
	 *
	 * @warning la ligne de l'élément recherché doit etre non négative
	 * et strictement inférieure aux nombre de lignes que contient le
	 * tableau.
	 *
	 * @warning la colone de l'élément recherché doit etre non négative
	 * et inférieure ou égale à la ligne de celui ci.
	 */
	inline Type & operator()(const int &  row , const int & column)
	{
		assert(row >= 0);
		assert(column >= 0);
		
		assert(row < getNbrOfRows());
		assert(column <= row);

		return getData(row,column);
	}
	
	/**
	 * Méthode pour changer la valeur d'un élément du tableau.
	 *
	 * @param value la valeur avec laquelle on souhaite remplir le containor.
	 *
	 * @param row la ligne de l'élément à modifier.
	 *
	 * @param column la colone de l'élément à modifier.
	 *
	 * @warning le numéro de la ligne doit etre positif ou nul et strictement 
	 * inférieur aux nombres de lignes que comporte le tableau.
	 * 
	 * @warning le numéro de la colone doit etre positif ou nul et inférieur ou 
	 * égale à la ligne.
	 */
	inline void setData(const int & row , const int & column , const Type & value)
	{
		assert(row >= 0);
		assert(row < getNbrOfRows());

		assert(column >= 0);
		assert(column <= row);

		containor<Type>::setData( hash(row,column) , value);
	}

	/**
	 * Méthode pour fixer les valeurs du tableau à partir d'un autre tableau.
	 *
	 * @param array le tableau que l'on souhaite recopier dans l'objet.
	 */
	inline void setDatas(const containorTri2D<Type> & array)
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
	 * Méthode pour modifier la taille du tableau.
	 *
	 * @param nbrOfRows le nombre de lignes (et de colones que pourra
	 * contenir le tableau.
	 *
	 * @warning il faut que le nombre de lignes soit strictement positif.
	 *
	 */
	inline void setSizes(const int & nbrOfRows)
	{
		assert(nbrOfRows > 0);
		int size =  (nbrOfRows * (nbrOfRows+1)) / 2;
		NbrOfRows = nbrOfRows;
		containor<Type>::setSizes(size);
	}

/**
 * Méthode pour afficher joliment les containor triangulaires bidimensionnels.
 * 
 * @param outStream le flux dans lequel on souhaite afficher le containor.
 *
 * @return void.	
 */
inline virtual void writePretty(ostream & outStream) const
{
	int row , column , nbrOfRows;

	nbrOfRows = getNbrOfRows();

	for( row = 0 ; row < nbrOfRows ; row++) {
		for( column = 0 ; column < nbrOfRows ; column++) {

			if(column <= row ) {
					outStream << setw(8) << getData(row , column);
				} else {
					outStream << setw(8) << "x"; 
				}

		}

		if(row < (nbrOfRows - 1) ) {
			outStream  << endl;
		}
		
	}

	
	return;
} 

};



#endif
