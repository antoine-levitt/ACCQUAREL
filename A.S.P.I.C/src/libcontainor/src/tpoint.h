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
#ifndef _T_POINT_
#define _T_POINT_

#include <iostream>
#include <string.h>
#include <assert.h>

using namespace std;

/**
 *  Classe pour manipuler des tableaux de taille
 * fixe (connue à la compilation).
 */
template<class Type , int Dimension>
class tpoint
{
private:

	/**
	 * Tableau de'elements de type Type est de taille Dimension.
	 */
	Type Datas[Dimension];

protected:

	/**
	 * Méthode pour copier le contenu d'un objet point
	 * dans l'objet courrant.
	 * 
	 * @param p le point que l'on souhaite copier.
	 * @return void.
	 */
	inline void copy(const tpoint<Type,Dimension> & p)
	{
		for(int i=0 ; i < Dimension ; i++) setData(i ,  p.getData(i)); 
	}

public:

	/**
	 * Constructeur par défaut.
	 */
	inline tpoint(void)
	{
		;
	}

	/**
	 * Constructeur avec la meme valeur pour tout le monde.
	 *
	 * @param value la valeur que prendrons tout les éléments du tableau
	 * lorsqu'il aura été créé.
	 */
	inline tpoint(const Type &  value)
	{
		setDatas(value);
	}

	/**
	 * Constrcuteur avec les valeur des trois éléments du tableau.
	 *
	 * @param x la valeur du premier element.
	 * @param y la valeur du second element.
	 * @param z la valeur du troisième élément. 
	 *
	 * @warning ce constructeur ne doit etre utilisé que pour 
	 * les objets de type tpoint<Type,3>.
	 */
	inline tpoint(const Type &  x , const Type & y , const Type & z)
	{
		assert( Dimension == 3);

		setDatas(x,y,z);
	}

	/**
	 * Constrcuteur avec les valeur des deux éléments du tableau.
	 *
	 * @param x la valeur du premier element.
	 * @param y la valeur du second element.
	 *
	 * @warning ce constructeur ne doit etre utilisé que pour 
	 * les objets de type tpoint<Type,2>.
	 */
	inline tpoint(const Type &  x , const Type &  y)
	{
		assert(Dimension == 2);

		setDatas(x,y);
	}

	/**
	 * Constructeur de copie.
	 *
	 * @param p le tableau à copier dans l'objet contruit.
	 */
	inline tpoint(const tpoint<Type,Dimension> & p)
	{
		copy(p);
	}

	/**
	 * Destructeur.
	 */
	inline virtual ~tpoint(void)
	{
		;
	}

	/**
	 * Méthode constante qui permet d'accéder un élément du tableau.
	 * 
	 * @param i le numéro de l'élément.
	 *
	 * @return une référence constante vers la i-ème coordonnée.
	 *
	 * @warning l'appel de cette fonction avec i négatif
	 * ou i supérieur ou égal à la dimension provoque l'arret 
	 * du programme.
	 */
	inline const Type & getData(int i) const
	{
		assert( i >= 0);
		assert( i < getDimension());

		return Datas[i];
	}

	/**
	 * Méthode qui permet d'accéder un élément du tableau.
	 * 
	 * @param i le numéro de l'élément.
	 *
	 * @return une référence vers le i-ème élément.
	 *
	 * @warning l'appel de cette fonction avec i négatif
	 * ou i supérieur ou égal à la dimension provoque l'arret 
	 * du programme.
	 */
	inline Type & getData(int i)
	{
		assert( i >= 0);
		assert( i < Dimension);

		return Datas[i];
	}

	/**
	 * Méthode pour connaitre la diùmension du tableau.
	 *
	 * @return la valeur de la dimension du tableau.
	 */
	inline const int getDimension(void) const
	{
		return Dimension;
	}

	/**
	 * Opérateur constant d'accès aux éléments du tableau.
	 *
	 * @param i le position de l'élément recherché.
   * 
	 * @return la une référence constante vers l'élément recherché.
	 *
	 * @warning la position de m'élément demandé doit etre positive ou
	 * nulle et strictement inférieure à la dimension du tableau.
	 */
	inline const Type & operator[] (const int & i) const
	{
		assert(i >= 0);
		assert(i < Dimension);

		return getData(i);
	}

	/**
	 * Opérateur constant d'accès aux éléments du tableau.
	 *
	 * @param i le position de l'élément recherché.
   * 
	 * @return la une référence constante vers l'élément recherché.
	 *
	 * @warning la position de m'élément demandé doit etre positive ou
	 * nulle et strictement inférieure à la dimension du tableau.
	 */
	inline Type & operator[] (const int & i)
	{
		assert(i >= 0);
		assert(i < Dimension);

		return getData(i);
	}

	/**
	 * Opérateur d'affectation.
	 *
	 * @param p le tableau que l'on souhaiote copier dans l'objet appelant.
	 *
	 * @return une référence vres le tableau modifié.
	 */
	inline tpoint<Type,Dimension> & operator= (const tpoint<Type,Dimension> & p)
	{
		copy(p);
		return *this;
	}

	/**
	 * Méthode pour caster l'objet tpoint<Type,Dimension> vers un pointeur de Type.
	 *
	 * @return un pointeur constant vers le tableau de données.
	 */
	inline operator const Type * (void) const 
	{
		return Datas;
	}

	/**
	 * Méthode qui permet la valeur d'un élément du tableau.
	 *
	 * @param i la position de l'élément dont on souhaite modifier la valeur.
	 * @param value la valeur à donner au i-ème élément du tableau.
	 *
	 * @warning l'appel de cette fonction avec i négatif
	 * ou i supérieur ou égal à la dimension provoque l'arret 
	 * du programme.
	 */
	inline void setData(const int & i , const Type & value)
	{
		assert( i >= 0);
		assert( i < Dimension);

		Datas[i] = value;
	}

	/**
	 * Méthode qui permet de donner à tous les éléments du 
	 * tableau une unique valeur.
	 *
	 * @param value la valeur à donner à tous les éléments du tableau.
	 */
	inline void setDatas(const Type & value)
	{
		for(int i=0 ; i < Dimension ; i++)
			setData(i,value);
	}

	/**
	 * Méthode pour donner des valeurs à tout les éléments
	 * d'un tableau de dimension 3 d'un seul coup.
	 *
	 * @param x la valeur du premier element.
	 * @param y la valeur du second element.
	 * @param z la valeur du troisième élément. 
	 *
	 * @warning cette méthode ne doit etre utilisé que pour 
	 * les objets de type tpoint<Type,3>.
	 */
	inline void setDatas(const Type & x , const Type & y , const Type & z)
	{
		assert(Dimension == 3);

		setData(0,x);
		setData(1,y);
		setData(2,z);
	}

	/**
	 * Méthode pour donner des valeurs à tout les éléments
	 * d'un tableau de dimension 3 d'un seul coup.
	 *
	 * @param x la valeur du premier element.
	 * @param y la valeur du second element.
	 *
	 * @warning cette méthode ne doit etre utilisé que pour 
	 * les objets de type tpoint<Type,2>.
	 */
	inline void setDatas(const Type & x , const Type & y) 
	{
		assert(Dimension == 2);

		setData(0,x);
		setData(1,y);
	}

	/**
	 * Méthode pour donner des valeurs à tout les éléments
	 * d'un tableau en les recopiant depuis un autre tableau.
	 *
	 * @param p le tableau à recopier.
	 */
	inline void setDatas(const tpoint<Type,Dimension> & p)
	{
		copy(p);
	}

	/**
	 * Méthode pour écrire le tableau dans un flux.
	 *
	 * @param out le flux dans lequel on souhaite écrire le tableau.
	 *
	 * @see operator<<(ostream & , const tpoint<Type,Dimension> & )
	 */
	inline virtual void write(ostream & out) const
	{
		out << "(";
		for(int i=0 ; i < Dimension  ; i++)  {
			out << getData(i) << (i < (Dimension-1) ? " , " : "" ) ;
		}
		out << ")";
	}

};	

/**
 * Opérateur écrire un tableau de type tpoint<Type,Dimension> dans 
 * flux.
 * 
 * @param out le flux dans lequel on souhaite écrire le tableau.
 *
 * @param p le tableau à écrire dans le flux.
 *
 * @return une référence vers le flux contenant la description du tableau.
 */
template<class Type , int Dimension>
inline ostream & operator<<(ostream & out , const tpoint<Type,Dimension> & p) 
{
	p.write(out);
	return out;
}

#endif



