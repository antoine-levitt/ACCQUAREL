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
#ifndef _CONTAINOR_
#define _CONTAINOR_

#include <assert.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

/**
 * @class containor.
 * Classe pour la gestion d'un conteneur d'un certains nombres d'objet.
 * 
 * Le propos de cette classe est de ne plus utiliser de new et delete
 * dans les programme, mais des objets containor.
 */
template<class Type>
class containor
{
	private:

	/**
	 * Le nombre de données qui sont presentes dans le conteneur.
	 */
	int NbrOfRows;
	
	/**
	 * Tableau de données.
	 */
	Type * Datas;

protected:

	/**
	 * Méthode qui permet d'allouer l'espace mémoire.
	 *
	 * @param nbrOfRows le nombre d'éléments qui pourront etre 
	 * contenues dans l'objet.
	 */
	inline void alloc(int nbrOfRows)
	{
		assert(NbrOfRows == 0);
		assert(nbrOfRows > 0);
		NbrOfRows = nbrOfRows;
		Datas = new Type[nbrOfRows];
	}

	/**
	 * Méthode qui permet de copier un conteneur dans l'objet courrant.
	 *
	 * @param c le containor à copier dans l'objet courrant.
	 */
	inline void copy(const containor<Type> & c)
	{
		int nbrOfRows = c.getNbrOfRows();
		setSizes(nbrOfRows);
		for(int i=0 ; i < nbrOfRows ; i++) setData(i , c.getData(i));
	}

	/**
	 * Méthode qui libère l'espace mémoire.
	 */
	inline void free(void)
	{
		if(NbrOfRows == 0) return;
		delete [] Datas;
		Datas = NULL;
		NbrOfRows = 0;
	}

public:

	/**
	 * Constructeur par défaut.
	 * 
	 * Construit un objet vide qui ne peut contenir aucune 
	 * donnée.
	 */
	inline containor(void)
		: NbrOfRows(0) , 
			Datas(NULL)
	{
		;
	}

	/**
	 * Constructeur avec une spécification de la taille.
	 *
	 * @param nbrOfRows la taille du conteneur qui doit etre construit.
	 *
	 * @warning il faut que le nombre d'élements du tableu soit strictement
	 * positif.
	 */
	inline containor(const int & nbrOfRows)
		: NbrOfRows(0) , 
			Datas(NULL)
	{
		assert (nbrOfRows > 0);

		this->setSize(nbrOfRows);
	}

	/**
	 * Constructeur de copie.
	 *
	 * @param c le containor à copier dans l'objet construit.
	 */
	inline containor(const containor<Type> & c) 
		: NbrOfRows(0) , 
			Datas(NULL)
	{
		copy(c);
	}

	/**
	 * Destrucuteur.
	 */
	inline virtual ~containor(void)
	{
		free();
	}

	/**
	 * Méthode pour remettre tout le conteneur à zéro.
	 *
	 * @warning ici remettre à zéro revient à créer un containor vide.
	 */
	inline void clear(void) 
	{
		free();
	}

	/**
	 * Methode pour savoir si l'objet que l'on manipule est vide.
	 *
	 * @return true lorsque l'objet est vide, false sinon.
	 */
	inline bool empty(void) const
	{
		if(getNbrOfRows() == 0) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Méthode pour connaitre la i-ème valeur contenue dans le tableau.
	 *
	 * @param i la position de la données à laquelle on souhaite accéder.
	 * @return une référence constante vers le i-ème élément du containor.
	 *
	 * @warning lorsque l'on passe un argument négatif ou supérieur
	 * au nombre de données cela provoque l'arret du programe.
	 */

	inline const Type & getData(const int & i) const
	{
		assert(i >= 0);
		assert(i < getNbrOfRows());

		return Datas[i];
	}

	/**
	 * Méthode pour connaitre la i-ème valeur contenue dans le tableau.
	 *
	 * @param i la position de la données à laquelle on souhaite accéder.
	 * @return une référence vers le i-ème élément du containor.
	 *
	 * @warning lorsque l'on passe un argument négatif ou supérieur
	 * au nombre de données cela provoque l'arret du programe.
	 */

	inline Type & getData(const int & i)
	{
		assert(i >= 0);
		assert(i < getNbrOfRows());
		return Datas[i];
	}

	/**
	 * Méthode accéder aux tableau comme si c'était un pointeur.
	 *
	 * @return le pointeur qui contient tout les éléments.
	 */
	inline const Type * getDatas(void) const
	{
		return Datas;
	}

	/**
	 * Méthode pour connaitre la taille du tableau.
	 * 
	 * @return le nombre d'éléments qui peuvent etre stockés
	 * dans l'objet.
   */
	inline const int & getNbrOfRows(void) const
	{
		return NbrOfRows;
	}

	/**
	 * Méthode pour connaitre la taille du tableau dans une direction donnée.
	 *
	 * Cette méthode est inutile, dans le cadre d'un containor tout simple, mais
	 * par homogénéité il me semble utile de la mettre. (Par ex. pour une future
	 * implémentation templatisée de ces classes).
	 * 
	 * @param dimension la dimension dans laquelle on souhaite connaitre la taille, 
	 * ici ce paramêtre doit etre 0.
	 *
	 * @return le nombre de lignes du le tableau.
	 *
	 * @warning la dimension doit prendre la valeur 0 ( pour accéder aux nombre de lignes).
	 */
	inline const int & getSize(const int & dimension) const
	{
		assert(dimension == 0);
		return getNbrOfRows();
	}

	/**
	 * Méthode pour connaitre la taille du tableau.
	 * 
	 * @return le nombre d'éléments qui peuvent etre stockés
	 * dans l'objet.
   */
	inline const int & getSizes(void) const
	{
		return getNbrOfRows();
	}

	/**
	 * Méthode pour connaitre la taille de l'espace mémoire alloué.
	 *
	 * @return le nombre d'éléments qui on étés alloués âr la méthode
	 * alloc.
	 */
	inline const int & getStorageSize(void) const
	{
		return getNbrOfRows();
	}

	/**
	 * Opérateur pour accéder aux éléments du conteneur.
	 *
	 * @param i la position de la données à laquelle on souhaite accéder.
	 * @return une référence constante vers le i-ème élément du containor.
	 *
	 * @warning lorsque l'on passe un argument négatif ou supérieur
	 * au nombre de données cela provoque l'arret du programe.
	 */
	inline const Type & operator[] (const int & i) const
	{
		assert(i >= 0);
		assert(i < getNbrOfRows());
		return getData(i);
	}

	/**
	 * Opérateur pour accéder aux éléments du conteneur.
	 *
	 * @param i la position de la données à laquelle on souhaite accéder.
	 * @return une référence vers le i-ème élément du containor.
	 *
	 * @warning lorsque l'on passe un argument négatif ou supérieur
	 * au nombre de données cela provoque l'arret du programe.
	 */
	inline Type & operator[] (const int & i) 
	{
		assert(i >= 0);
		assert(i < getNbrOfRows());
		return getData(i);
	}

	/**
	 * Opérateur de copie
	 *
	 * @param c le containor que l'on souhaite copier dans l'objet courrant.
	 * @return une référence vers le nouvel objet qui contient le containor c.
	 */

	inline containor & operator=(const containor<Type> & c)
	{
		copy(c);
		return *this;
	}	

	/**
	 * Méthode pour changer la valeur d'un élément du tableau.
	 *
	 * @param value la valeur avec laquelle on souhaite remplir le containor.
	 *
	 * @param i la positon de l'élément à modifier.
	 *
	 * @warning lorsque l'on passe comme position un argument négatif ou supérieur
	 * au nombre de données cela provoque l'arret du programe.
	 */
	inline void setData(const int & i , const Type & value)
	{
		assert(i >= 0);
		assert(i < getNbrOfRows());

		Datas[i] = value;
	}

	/**
	 * Méthode pour remplir le conteneur avec une seule valeure.
	 *
	 * @param value la valeur avec laquelle on souhaite remplir le containor.
	 */
	inline void setDatas(const Type & value)
	{
		for(int i=0 ; i < getNbrOfRows() ; i++) setData(i , value);
	}

	
	inline void setDatas(const containor<Type> & array)
	{
		copy(array);
	}


	/**
	 * Méthode pour spécifier la taille du containor.
	 *
	 * @param nbrOfRows le nombre de données qui pourront etre contenues
	 * dans l'objet.
	 *
	 * @warning le nombre d'éléments que peut contenir le tableau doit etre un
	 * nombre strictement positif.
	 *
	 * @see clear.
	 */
	inline void setSizes(const int & nbrOfRows)
	{
		assert (nbrOfRows > 0);

		if(nbrOfRows == getNbrOfRows()) return;
		if(getNbrOfRows() > 0) free();
		alloc(nbrOfRows);
	}

	/**
	 * Méthode qui permet d'afficher le conteneur sous une forme
	 * lisible pour l'homme.
	 * 
	 * Il est iportant de remarquer que cette méthode est virtuelle. 
	 * Cela permet à tous les objets qui héritent de la classe containor
	 *  et qui  réimplémentent cette méthode, utiliser
	 * l'operator<< qui n'est qu'une encapsulation de cette méthode.
	 * 
	 * @param out le flux dans lequel on veut écrire le contrainor. 
	 */
	inline virtual void writePretty(ostream & out) const
	{
		int i , nbrOfRows = getNbrOfRows();
		out << "[ " ;
		for( i = 0 ; i < nbrOfRows ; i++) {
			out << getData(i) << ( i == (nbrOfRows -1) ? " ]" : " , " );
		}
	}

	/**
	 * Méthode qui permet d'afficher le conteneur sous une forme
	 * lisible pour l'homme.
	 * 
	 * @param out le flux dans lequel on veut écrire le contrainor.
	 * @param enclosingTagName le nom du tag qui marquera le début et
	 * la fin des données du containor.
	 *
	 */
	inline void writeXML(ostream & out , const string & enclosingTagName = "Containor") const
	{
		int i , nbrOfRows = getNbrOfRows();
		
		out << "<" << enclosingTagName << ">" << endl;
		
		for( i = 0 ; i < nbrOfRows ; i++) {
			out << "<Coefficient>" << endl;				
			out << "<Row>" << i << "</Row>" << endl;
			out << "<Value>" << getData(i)<< "</Value>" << endl; 
			out << "</Coefficient>" << endl;	
		}

		out << "</" << enclosingTagName << ">" << endl;
	}

};

/**
 * Opérateur d'écriture d'un containor dns un flux.
 *
 * @param out le flux dans lequel on veut écrire le containor.
 * @param c le containor à écrire dans le flux.
 *
 * @param out le flux dans lequel on souhaite écrire le containor.
 */
template<class Type>
inline ostream & operator << ( ostream & out , const containor<Type> & c) 
{
	c.writePretty(out);
	return out;
}

#endif
