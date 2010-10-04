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
#ifndef _containor3D_
#define _containor3D_

#include "containor.h"
#include <iostream>
#include <iomanip>
#include "ipoint.h"
using namespace std;


/**
 * La classe des tableaux tri dimensionnels.
 */
template<typename Type>
class containor3D : public containor<Type>
{
private:
	
	/**	 
	 * Taille de la première dimension.
	 */
	int NbrOfColumns;
	
	/**
	 * Taille de la seconde dimension.
	 */
	int NbrOfLayers;
	
	/**
	 * Taille de la troisème dimension.
	 */
	int NbrOfRows;
	
protected:
	
	/**
	 * Méthode qui permet de copier un objet de de type containor3D<Type> dans l'objet courrant.
	 *
	 * @param c le cube à copier.
	 */
	inline void copy(const containor3D<Type> & c)
	{

		if(c.empty()) {
			clear();
			return;
		}

		NbrOfRows = c.getNbrOfRows();
		NbrOfColumns = c.getNbrOfColumns();
		NbrOfLayers = c.getNbrOfLayers();
		
		containor<Type>::copy(c);
	}

	/**
	 * Méthode qui permet de retrouver la position d'un élément du tableau
	 * tri dimensionel dans le vectuer de stockage.
	 *
	 * @param row la ligne de l'élément.
	 *
	 * @param column la colone de l'élément.
	 *
	 * @param layer la couche de l'élément.
	 *
	 * @return la position de cet élément dans le tableau qui sert à stoquer les
	 * valeurs des éléments.
	 *
	 * @warning il faut que le paramètre row soit non négatif et strictement
	 * inférieur au nombre de lignes de l'objet.
	 *
	 * @warning il faut que le paramètre column soit non négatif et strictement
	 * inférieur au nombre de colones de l'objet.
	 *
	 * @warning il faut que le paramètre layer soit non négatif et strictement
	 * inférieur au nombre de couches de l'objet.
	 */
	inline int hash(const int & row , const int & column , const int & layer) const
	{
		assert(row >= 0);
		assert(row < getNbrOfRows());
		
		assert(column >= 0);
		assert(column < getNbrOfColumns());
		
		assert(layer >= 0);
		assert(layer < getNbrOfLayers());

		int h;
		
		h= row + getNbrOfRows() * ( column + getNbrOfColumns() * layer);
		
		return h;
	}


public:

	/**
	 * Constructeur par défaut.
	 *
	 * Ce constructeur construit un objet totalement vide.
	 */
	inline containor3D(void)
		: containor<Type>() ,
		NbrOfColumns(0) , 
		NbrOfLayers(0) ,
		NbrOfRows(0)
	{
		;
	}

	/**
	 * Constructeur avec spécification des tailles.
	 * 
	 * @param nbrOfRows le nombre de lignes que contiendra l'objet.
	 *
	 * @param nbrOfColumns le nombre de colones que contiendra l'objet.
	 *
	 * @param nbrOfLayers le nombre de couches que contiendra l'objet.
	 *
	 * @warning il faut que le paramètre nbrOfRows soit strictement positif.
	 *
	 * @warning il faut que le paramètre nbrOfColumns soit strictement positif.
	 *
	 * @warning il faut que le paramètre nbrOfLayer soit strictement positif.
	 */
	inline containor3D(const int & nbrOfRows , const int & nbrOfColumns , const int & nbrOfLayers)
		: containor<Type>() ,
		NbrOfColumns(0) , 
		NbrOfLayers(0) ,
		NbrOfRows(0)
	{
		this->setSizes(nbrOfRows , nbrOfColumns , nbrOfLayers);
	}

	/**
	 * Constructeur de copie.
	 *
	 * @param c le containor à copier d'ans l'objet nouvellement construit.
	 */
	inline containor3D(const containor3D<Type> & c)
		: containor<Type>() ,
		NbrOfColumns(0) , 
		NbrOfLayers(0) ,
		NbrOfRows(0)
	{
		copy(c);
	}

	/**
	 * Destructeur.
	 */
	virtual ~containor3D(void)
	{
		clear();
	}

	/**
	 * Méthode qui permet de nettoyer le tableau.
	 *
	 * @warning après l'appel de cette méthode, l'objet ne
	 * appelant ne peut plus contenir une seule donnée.
	 */
	inline void clear(void)
	{
		NbrOfRows = 0;
		NbrOfColumns = 0;
		NbrOfLayers = 0;

		containor<Type>::clear();
	}
	
	/**
	 * Méthode pour accéder à un élément du tableau.
	 *
	 * @param row la ligne de l'élément.
	 *
	 * @param column la colone de l'élément.
	 *
	 * @param layer la couche de l'élément.
	 *
	 * @return une reférence vers l'élément du tableau à la position (row , column , layer).
	 *
	 * @warning il faut que le paramètre row soit non négatif et strictement
	 * inférieur au nombre de lignes de l'objet.
	 *
	 * @warning il faut que le paramètre column soit non négatif et strictement
	 * inférieur au nombre de colones de l'objet.
	 *
	 * @warning il faut que le paramètre layer soit non négatif et strictement
	 * inférieur au nombre de couches de l'objet.
	 */
	inline Type & getData(const int & row , const int & column , const int & layer) 
	{
		assert(row >= 0);
		assert(row < getNbrOfRows());
		
		assert(column >= 0);
		assert(column < getNbrOfColumns());
		
		assert(layer >= 0);
		assert(layer < getNbrOfLayers());
		
		return containor<Type>::getData(hash(row , column , layer));
	}

	/**
	 * Méthode constante pour accéder à un élément du tableau.
	 *
	 * @param row la ligne de l'élément.
	 *
	 * @param column la colone de l'élément.
	 *
	 * @param layer la couche de l'élément.
	 *
	 * @return une référence cosntante vers l'élément du tableau à la position (row , column , layer).
	 *
	 * @warning il faut que le paramètre row soit non négatif et strictement
	 * inférieur au nombre de lignes de l'objet.
	 *
	 * @warning il faut que le paramètre column soit non négatif et strictement
	 * inférieur au nombre de colones de l'objet.
	 *
	 * @warning il faut que le paramètre layer soit non négatif et strictement
	 * inférieur au nombre de couches de l'objet.
	 */
	inline const Type & getData(const int & row , const int & column , const int & layer) const 
	{
		assert(row >= 0);
		assert(row < getNbrOfRows());
		
		assert(column >= 0);
		assert(column < getNbrOfColumns());
		
		assert(layer >= 0);
		assert(layer < getNbrOfLayers());
		
		return containor<Type>::getData(hash(row , column , layer));
	}

	/**
	 * Méthode pour accéder aux tailles du tableau dans les trois dimensions.
	 *
	 * @return un objet de type ipoint<3> dont la première composante contient le nombre
	 * de lignes , la seconde le nombre de colones et la dernière le nombres de couches.
	 */
	inline ipoint<3> getSizes(void) const
	{
		ipoint<3> sizes;

		sizes[0] = getNbrOfRows();
		sizes[1] = getNbrOfColumns();
		sizes[2] = getNbrOfLayers();

		return sizes;
	}
	
	/**
	 * Méthode pour accéder à la taille du tableau dans l'une des trois dimensions.
	 *
	 * @param dim la dimension dans laquelle on souhaite connaitre la taille. 
	 *
	 * @return Lorsque le paramètre dim prend la valeur 0, la méthode renvoie le nombre
	 * de lignes ; lorsque le paramètre prend la valeur 1, la méthode renvoie le nombre 
	 * de colones ; enfin pour une valeur de 2 en entrée alors la méthode renoie le nombre
	 * de couches.
	 *
	 * @warning le paramètre dim ne peut prendre que les valeurs 0 , 1 , ou 2.
	 */
	inline virtual const int & getSize(const int & dim) const
	{
		assert( dim >= 0);
		assert(dim < 3);

		if(dim==0) {
			return getNbrOfRows();
		} else if(dim==1) {
			return getNbrOfColumns();
		} else {
			return getNbrOfLayers();
		}
	}

	/**
	 * Méthode pour accéder aux nombres de lignes que contient le tableau.
	 *
	 * @return le nombre de lignes que contient le tableau.
	 */
	inline const int & getNbrOfRows(void) const
	{
		return NbrOfRows;
	}

	/**
	 * Méthode pour accéder aux nombres de colones que contient le tableau.
	 *
	 * @return le nombre de colones que contient le tableau.
	 */
	inline const int & getNbrOfColumns(void) const
	{
		return NbrOfColumns;
	}

	/**
	 * Méthode pour accéder aux nombres de couches que contient le tableau.
	 *
	 * @return le nombre de couches que contient le tableau.
	 */
	inline const int & getNbrOfLayers(void) const
	{
		return NbrOfLayers;
	}

	/**
	 * Opérateur d'affectation.
	 */
	inline containor3D<Type> & operator= (const containor3D<Type> & c)
	{
		copy(c);
		return *this;
	}

	/**
	 * Opérateur constant pour accéder à un élément du tableau.
	 *
	 * @param row la ligne de l'élément.
	 *
	 * @param column la colone de l'élément.
	 *
	 * @param layer la couche de l'élément.
	 *
	 * @return une référence cosntante vers l'élément du tableau à la position (row , column , layer).
	 *
	 * @warning il faut que le paramètre row soit non négatif et strictement
	 * inférieur au nombre de lignes de l'objet.
	 *
	 * @warning il faut que le paramètre column soit non négatif et strictement
	 * inférieur au nombre de colones de l'objet.
	 *
	 * @warning il faut que le paramètre layer soit non négatif et strictement
	 * inférieur au nombre de couches de l'objet.
	 */
	inline const Type & operator()(const int & row , const int & column , const int & layer) const
	{
		assert(row >= 0);
		assert(row < getNbrOfRows());
		
		assert(column >= 0);
		assert(column < getNbrOfColumns());
		
		assert(layer >= 0);
		assert(layer < getNbrOfLayers());

		return getData(row , column , layer);
	}

	/**
	 * Opérateur pour accéder à un élément du tableau.
	 *
	 * @param row la ligne de l'élément.
	 *
	 * @param column la colone de l'élément.
	 *
	 * @param layer la couche de l'élément.
	 *
	 * @return une référence vers l'élément du tableau à la position (row , column , layer).
	 *
	 * @warning il faut que le paramètre row soit non négatif et strictement
	 * inférieur au nombre de lignes de l'objet.
	 *
	 * @warning il faut que le paramètre column soit non négatif et strictement
	 * inférieur au nombre de colones de l'objet.
	 *
	 * @warning il faut que le paramètre layer soit non négatif et strictement
	 * inférieur au nombre de couches de l'objet.
	 */
	inline Type & operator()(const int & row , const int & column , const int & layer)
	{
		assert(row >= 0);
		assert(row < getNbrOfRows());
		
		assert(column >= 0);
		assert(column < getNbrOfColumns());
		
		assert(layer >= 0);
		assert(layer < getNbrOfLayers());

		return getData(row , column , layer);
	}

	/**
	 * Opérateur constant pour accéder à un élément du tableau.
	 *
	 * @param position le triplet d'entier qui contient la ligne , la colone et la couche
	 * de l'élément dont on souhaite conniatre la valeur.
	 *
	 * @return une référence cosntante vers l'élément du tableau à la position (position[0] , position[1] , position[2]).
	 *
	 * @warning il faut que le paramètre position[0] soit non négatif et strictement
	 * inférieur au nombre de lignes de l'objet.
	 *
	 * @warning il faut que le paramètre position[1] soit non négatif et strictement
	 * inférieur au nombre de colones de l'objet.
	 *
	 * @warning il faut que le paramètre position[2] soit non négatif et strictement
	 * inférieur au nombre de couches de l'objet.
	 */
	inline const Type & operator()(const ipoint<3> & position) const
	{
		assert(position[0] >= 0);
		assert(position[0] < getNbrOfRows());
		
		assert(position[1] >= 0);
		assert(position[1] < getNbrOfColumns());
		
		assert(position[2] >= 0);
		assert(position[2] < getNbrOfLayers);
		
		return getData(position[0], position[1], position[2]);
	}

	/**
	 * Opérateur pour accéder à un élément du tableau.
	 *
	 * @param position le triplet d'entier qui contient la ligne , la colone et la couche
	 * de l'élément dont on souhaite conniatre la valeur.
	 *
	 * @return une référence vers l'élément du tableau à la position (position[0] , position[1] , position[2]).
	 *
	 * @warning il faut que le paramètre position[0] soit non négatif et strictement
	 * inférieur au nombre de lignes de l'objet.
	 *
	 * @warning il faut que le paramètre position[1] soit non négatif et strictement
	 * inférieur au nombre de colones de l'objet.
	 *
	 * @warning il faut que le paramètre position[2] soit non négatif et strictement
	 * inférieur au nombre de couches de l'objet.
	 */
	inline Type & operator()(const ipoint<3> & position)
	{
		assert(position[0] >= 0);
		assert(position[0] < getNbrOfRows());
		
		assert(position[1] >= 0);
		assert(position[1] < getNbrOfColumns());
		
		assert(position[2] >= 0);
		assert(position[2] < getNbrOfLayers());
		
		return getData(position[0], position[1], position[2]);
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
	 * @param layer la couche de l'élément à modifier.
	 *
	 * @warning le numéro de la ligne doit etre positif ou nul et strictement 
	 * inférieur aux nombres de lignes que comporte le tableau.
	 * 
	 * @warning le numéro de la colone doit etre positif ou nul et inférieur
	 * strictement aux nombres de colones du tableau.
	 *
	 * @warning le numéro de la couche doit etre positif ou nul et inférieur
	 * strictement au nombre de couches du tableau.
	 */
	inline void setData(const int & row , const int & column , const int & layer , const Type & value)
	{
		assert(row >= 0);
		assert(row < getNbrOfRows());

		assert(column >= 0);
		assert(column <= getNbrOfColumns());

		assert(layer >= 0);
		assert(layer <= getNbrOfLayers());

		containor<Type>::setData(hash(row,column,layer) , value);
	}

	/**
	 * Méthode pour fixer les valeurs du tableau à partir d'un autre tableau.
	 *
	 * @param array le tableau que l'on souhaite recopier dans l'objet.
	 */
	inline void setDatas(const containor3D<Type> & array)
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
	 * Méthode pour fixer la taille du tableau.
	 *
	 * @param nbrOfRows le nombre de lignes du tableau.
	 * @param nbrOfColumns le nombre de colones du tableau.
	 * @param nbrOfLayers le nombre de coiuches du tableau.
   *
	 * @warning le nombre de lignes spécifié en argument doit etre strictement positif.
	 * @warning le nombre de colones spécifié en argument doit etre strictement positif.
	 * @warning le nombre de couches spécifié en argument doit etre strictement positif.
	 */
	inline void setSizes(const int & nbrOfRows , const int & nbrOfColumns , const int & nbrOfLayers)
	{
		assert(nbrOfRows > 0);
		assert(nbrOfColumns > 0);
		assert(nbrOfLayers > 0);

		NbrOfRows = nbrOfRows;
		NbrOfColumns = nbrOfColumns;
		NbrOfLayers = nbrOfLayers;
		
		containor<Type>::setSizes(nbrOfRows * nbrOfColumns * nbrOfLayers);
	}

	/**
	 * Méthode pour fixer la taille du tableau.
	 *
	 * @param sizes un triplet d'entier qui contient le nombre de lignes , le nombre de colones et le nombre de couches.
	 *
	 * @warning le nombre de lignes spécifié en argument doit etre strictement positif.
	 * @warning le nombre de colones spécifié en argument doit etre strictement positif.
	 * @warning le nombre de couches spécifié en argument doit etre strictement positif.
	 */
	inline void setSizes(const ipoint<3> & sizes)
	{
		assert(sizes[0] >= 0);
		assert(sizes[1] >= 0);
		assert(sizes[2] >= 0);
		
		setSizes(sizes[0] , sizes[1] , sizes[2]);
	}

	/**
	 * Méthode pour afficher joliment les containor 3D.
	 * 
   * @param outStream le flux dans lequel on souhaite afficher le containor.
   *
   * @return void.	
   */
	inline virtual void writePretty(ostream & outStream) const
	{
		int row , column , layer , nbrOfRows , nbrOfColumns , nbrOfLayers;

		nbrOfLayers = getNbrOfLayers();
		nbrOfRows = getNbrOfRows();
		nbrOfColumns = getNbrOfColumns();

		for(layer = 0 ; layer < nbrOfLayers ; layer ++) {
			outStream << setw(8*nbrOfColumns) << setfill('-') << "" << setfill(' ') << endl;
	
			
				for( row = 0 ; row < nbrOfRows ; row++) {
				for( column = 0 ; column < nbrOfColumns ; column++) {
					
					outStream << setw(8) << getData(row,column,layer);

				}

				if(row < (nbrOfRows - 1) ) {
					outStream  << endl;
				}

			}

			if(layer < (nbrOfLayers - 1) ) {
				outStream  << endl;
			}
		}
		return;
	} 
};

#endif

