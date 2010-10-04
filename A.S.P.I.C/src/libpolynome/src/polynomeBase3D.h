/* 
* The polynome library of the A.S.P.I.C. 
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
#ifndef _POLYNOME_BASE_3D_
#define _POLYNOME_BASE_3D_

#define _USES_MATH_DEFINES
#include <containorSparse3D.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include <string>
#include "string4Polynomes.h"
using namespace std;

/**
 * Classe qui permet de gérer les polynomes 
 * dans le cas tri - dimensionnel.
 */
template<typename Type>
class polynomeBase3D : public containorSparse3D<Type>
{
private:

	
protected : 
	
	/**
	 * Méthode pour afficher une fonction de base
	 *
	 * Le problème est de pouvoir afficher tout les polynomes
	 * sans avoir à réécrire de méthode write ou d'opérateur <<.
	 * Pour cela on demande à toutes les classes enfants de réimplémenter
	 * cette méthode comme cela chaque classe affichera correctement
	 * la base de polynômes qui lui convient.
	 *
	 * Par exemple dans polynome_3D le getBaseString affichera quelque
	 * chose comme X^dx . Y^dy . Z^dz alors que pour un polynome centré
	 * on pourra afcher (X-CenterX)^dx . (Y-CenterY)^dy . (Z-CenterZ)^dz
	 *
	 * @param degree le degréé de la fonction de base du polynome dont 
	 * on souhaite connaitre la représentation.
	 */
	virtual const string getBaseString(const ipoint<3> & degree) const;

		/**
	 * Méthode de base pour faire l'addition de deux polynomes.
	 *
	 * Cette méthode permet de calculer la somme des polynomes de base
	 * p et q et de stoquer le résultat de l'opération.
	 *
	 * @param p le premier polynome.
	 *
	 * @param q le second polynome.
	 *
	 * @warning cette méthode est brute de fonderie, c'est a dire que les
	 * classes enfants auront la responsabilité de la rendre accesible à 
	 * l'utilisateur.
	 */
	void add(const polynomeBase3D<Type> & p , const polynomeBase3D<Type> & q);

	/**
	 * Méthode de base pour faire l'addition unaire de deux polynomes.
	 *
	 * Cette méthode permet de calculer la somme des polynomes de base
	 * q et du polynome appelant de stoquer le résultat de l'opération dans ce
	 * meme polynome.
	 *
	 * @param q à ajouter.
	 *
	 * @warning cette méthode est brute de fonderie, c'est a dire que les
	 * classes enfants auront la responsabilité de la rendre accesible à 
	 * l'utilisateur.
	 */
	void add(const polynomeBase3D<Type> & q);

	/**
	 * Méthode pour faire la division d'un polynome par un scalaire.
	 *
	 * Cette méthode permet de calculer la division d'un polynome de base par un scalaire.
	 * Le résultat est stocké dans l'objet.
	 *
	 * @param poly le polynome dont on souhaite diviser tout les coefficients.
	 * 
	 * @param scalar le salaire.
	 *
	 * @warning cette méthode est brute de fonderie, c'est a dire que les
	 * classes enfants auront la responsabilité de la rendre accesible à 
	 * l'utilisateur.
	 */
	void divide(const polynomeBase3D<Type> & poly , const Type & scalar);

	/**
	 * Méthode pour faire la division unaire d'un polynome et d'un scalaire.
	 *
	 * Cette méthode permet de calculer la division de l'objet par un scalaire.
	 * Le résultat est stocké dans l'objet.
	 *
	 * @param scalar le salaire.
	 *
	 * @warning cette méthode est brute de fonderie, c'est a dire que les
	 * classes enfants auront la responsabilité de la rendre accesible à 
	 * l'utilisateur.
	 *
	 * @warning le scalaire doit etre non nul.
	 */
	void divide(const Type & scalar);

	/**
	 * Méthode pour faire la multiplication de deux polynomes.
	 *
	 * Cette méthode permet de calculer le produit des polynomes de base
	 * p et q et de stoquer le résultat de l'opération.
	 *
	 * @param p le premier polynome.
	 *
	 * @param q le second polynome.
	 *
	 * @warning cette méthode est brute de fonderie, c'est a dire que les
	 * classes enfants auront la responsabilité de la rendre accesible à 
	 * l'utilisateur.
	 */
	void multiply(const polynomeBase3D<Type> & p , const polynomeBase3D<Type> & q);
	
	/**
	 * Méthode pour faire la multiplication de deux polynomes.
	 *
	 * Cette méthode permet de calculer le produit des polynomes de base
	 * p et q et de stoquer le résultat de l'opération.
	 *
	 * @param q le second polynome.
	 *
	 * @warning cette méthode est brute de fonderie, c'est a dire que les
	 * classes enfants auront la responsabilité de la rendre accesible à 
	 * l'utilisateur.
	 */
	void multiply(const polynomeBase3D<Type> & q);

	/**
	 * Méthode pour faire la multiplication d'un polynome et d'un scalaire.
	 *
	 * Cette méthode permet de calculer le produit de l'objet avec un scalaire.
	 * Le résultat est stocké dans l'objet.
	 *
	 * @param poly le polynome dont on souhaite multiplier tout les coefficients
	 ù par un scalaire.
	 *
	 * @param scalar le salaire.
	 *
	 * @warning cette méthode est brute de fonderie, c'est a dire que les
	 * classes enfants auront la responsabilité de la rendre accesible à 
	 * l'utilisateur.
	 */
	void multiply(const polynomeBase3D<Type> & poly , const Type & scalar);

	/**
	 * Méthode pour faire la multiplication unaire d'un polynome et d'un scalaire.
	 *
	 * Cette méthode permet de calculer le produit de l'objet avec un scalaire.
	 * Le résultat est stocké dans l'objet.
	 *
	 * @param scalar le salaire.
	 *
	 * @warning cette méthode est brute de fonderie, c'est a dire que les
	 * classes enfants auront la responsabilité de la rendre accesible à 
	 * l'utilisateur.
	 */
	void multiply(const Type & scalar);

	/**
	 * Méthode de base pour faire l'addition de deux polynomes.
	 *
	 * Cette méthode permet de calculer la somme des polynomes de base
	 * p et q et de stoquer le résultat de l'opération dans l'objet qui 
	 * effectue l'opération.
	 *
	 * @param p le premier polynome.
	 *
	 * @param q le second polynome.
	 *
	 * @warning cette méthode est brute de fonderie, c'est a dire que les
	 * classes enfants auront la responsabilité de la rendre accesible à 
	 * l'utilisateur.
	 */
	void soustract(const polynomeBase3D<Type> & p , const polynomeBase3D<Type> & q);

	/**
	 * Méthode de base pour faire la soustraction unaire de deux polynomes.
	 *
	 * Cette méthode permet de calculer la différence du polynome et du polynome passé en argument.
	 * Le résultat de l'opération dans l'obet qui effectue l'opération.
	 *
	 * @param q à ajouter.
	 *
	 * @warning cette méthode est brute de fonderie, c'est a dire que les
	 * classes enfants auront la responsabilité de la rendre accesible à 
	 * l'utilisateur.
	 */
	void soustract(const polynomeBase3D<Type> & q);

public:
	
	/**
	 * Constructeur par défaut.
	 */
	polynomeBase3D(void);
	
	/**
	 * Constrcuteur de copie.
	 *
	 * @param poly le polynome avec lequel on souhaite initialiser l'objet construit.
	 */
	polynomeBase3D(const polynomeBase3D<Type> & poly);

	/**
	 * Destrcuteur.
	 */
	virtual ~polynomeBase3D(void);

	/**
	 * Méthode pour connaitre les plus grands degrés d'un polynome pour chacune des
	 * trois composantes.
	 *
	 * Cette méthode renvoie un triplet d'entiers dont la première composante 
	 * est le plus haut degré pour la variable x, la seconde le plus haut degré pour 
	 * la variable y et la dernière le plus haut degré pour la variable z. 
	 * 
	 * Par exemple pour le polynome x.y.z cette méthode renvoie (1,1,1) et non pas 3. 
	 *
	 * @warning cette methode ne renvoie pas le degré du polynome.
	 */
	int getDegree(void) const;

	/**
	 * Méthode pour connaitre les plus grands degrés d'un polynome pour chacune des
	 * trois composantes.
	 *
	 * Cette méthode renvoie un triplet d'entiers dont la première composante 
	 * est le plus haut degré pour la variable x, la seconde le plus haut degré pour 
	 * la variable y et la dernière le plus haut degré pour la variable z. 
	 * 
	 * Par exemple pour le polynome x.y.z cette méthode renvoie (1,1,1) et non pas 3. 
	 * @see getExtermPosition(void)
	 * @warning cette methode ne renvoie pas le degré du polynome.
	 */
	ipoint<3> getUniformMax4Degree(void) const;

	/**
	 * Opérateur d'assertion.
	 *
	 * @param poly le polynome à copier dans l'objet.
	 */
	polynomeBase3D & operator = ( const polynomeBase3D<Type> & poly);

	/**
	 * Méthode pour écrire le polynome dans un flux.
	 *
	 * @param os le lfux danslequel on souahite écrire le polynome.
	 */
	void write(ostream & os) const;
};

/**
 * Opérateur pour écrire un polynome_base_3D dans un flux.
 *
 * @param out le flux dans lequel on souhaite écrire le polynôme.
 *
 * @param poly le polynome à écrire dans le flux.
 *
 * @return le flux contenant le polynome.
 */
template<typename Type>
extern ostream & operator<<(ostream & out , const polynomeBase3D<Type> & poly);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Le constructeur par défaut.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
polynomeBase3D<Type>::polynomeBase3D(void)
: containorSparse3D<Type>()
{
  this->clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Le constructeur copie.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
polynomeBase3D<Type>::polynomeBase3D(const polynomeBase3D<Type> & poly)
: containorSparse3D<Type>()
{
  containorSparse3D<Type>::copy(poly);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Le destructeur.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
polynomeBase3D<Type>::~polynomeBase3D(void)
{
  this->clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire l'addition de deux polynomes.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase3D<Type>::add(const polynomeBase3D<Type> & p , const polynomeBase3D<Type> & q)
{
	containorSparse3D<Type>::copy(p);
	add(q);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire l'addition unaire de deux polynomes.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase3D<Type>::add(const polynomeBase3D<Type> & q)
{
	ipoint<3> degree;
	Type value;

	for(degree = q.begin() ; degree != q.end() ; degree = q.next(degree)) {
		value = q(degree) + this->getData(degree);
		setData(degree , value);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la division d'un polynome et d'un scalaire.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase3D<Type>::divide(const polynomeBase3D<Type> & poly , const Type & scalar)
{
	assert(scalar != 0);

	containorSparse3D<Type>::copy(poly);
	divide(scalar);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la division unaire d'un polynome et d'un scalaire.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase3D<Type>::divide(const Type & scalar)
{
	assert(scalar != 0);

	ipoint<3> i;
	Type value;
	
	for(i=this->begin() ; i != this->end() ; i= this->next(i)) {
		value = this->getData(i) / scalar;
		setData(i,value);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour conaitre le degré du polynome.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
int polynomeBase3D<Type>::getDegree(void) const
{
	/* TODO */
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour conaitre le degré du polynome.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
ipoint<3> polynomeBase3D<Type>::getUniformMax4Degree(void) const
{
	return this->getExtremPositions();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour conaitre le degré du polynome.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
const string polynomeBase3D<Type>::getBaseString(const ipoint<3> & degree) const
{
	return monome3D2string(degree);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Methode pour faire la multiplication de deux polynomes
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase3D<Type>::multiply(const polynomeBase3D<Type> & p , const polynomeBase3D<Type> & q)
{
	ipoint<3> i,j;
	Type value;

	this->clear();

	for(i=q.begin() ; i != q.end() ; i=q.next(i)) {
		for(j=p.begin() ; j != p.end() ; j=p.next(j)) {
			value = this->getData(i+j) + q(i) * p(j);
			setData(i+j , value);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Methode pour faire la multiplication unaire de deux polynomes
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase3D<Type>::multiply(const polynomeBase3D<Type> & q)
{
	polynomeBase3D tmp(*this);
	multiply(tmp,q);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la multiplication unaire de d'un polynome et d'un scalaire.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase3D<Type>::multiply(const Type & scalar)
{
	ipoint<3> degree;
	Type value;

	if(scalar == 0) {
		this->clear();
		return;
	}
	
	for(degree = this->begin() ; degree != this->end() ; degree = this->next(degree)) {
		value = scalar *  this->getData(degree);
		setData(degree , value);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la multiplication de d'un polynome et d'un scalaire.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
 void polynomeBase3D<Type>::multiply(const polynomeBase3D<Type> & p , const Type & scalar)
{
	containorSparse3D<Type>::copy(p);
	multiply(scalar);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la soustraction de deux polynomes.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase3D<Type>::soustract(const polynomeBase3D<Type> & p , const polynomeBase3D<Type> & q)
{
	containorSparse3D<Type>::copy(p);
	soustract(q);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la soustraction unaire de deux polynomes.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase3D<Type>::soustract(const polynomeBase3D<Type> & q)
{
	ipoint<3> degree;
	Type value;

	for(degree = q.begin() ; degree != q.end() ; degree = q.next(degree)) {
		value = this->getData(degree) - q(degree);
		setData(degree , value);
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method WRITE for the polynome.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase3D<Type>::write(ostream & os) const
{
	bool has_previous = false;
	Type coefficient;
	ipoint<3> i;

	if(this->empty()) {
		os << "0.";
		return;
	}

	for( i=this->begin() ; i != this->end() ; i = this->next(i)) {

		coefficient = this->getData(i);

		if(has_previous) {
			os << (coefficient > 0 ? " +" :" ");
		}

		has_previous = false;

		if(i==ipoint<3>(0)) {
			os << coefficient;
			has_previous = true;

		} else { 

			if(fabs(coefficient) != 1) {
				os << coefficient;
				has_previous = true;		
			} else if(coefficient == -1) {
				os << "-" ;
			}

			if(has_previous)
				os << "*";

			os  << getBaseString(i);
			has_previous = true;

		} 
	} 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Extern Out Stream operator.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
ostream & operator << (ostream & out , const polynomeBase3D<Type> & poly)
{
	poly.write(out);
	return out;
}
#endif

