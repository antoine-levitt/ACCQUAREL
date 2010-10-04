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
#ifndef _POLYNOME_BASE_
#define _POLYNOME_BASE_

#include <assert.h>
#include <containorSparse.h>
#include <iostream>
#include <map>
#include <string>
#include "string4Polynomes.h"

using namespace std;

/**
 * Cette classe permet de stocker les coefficients
 * d'un objet sur une base qui peut etre 
 * inifinie.
 */
template<typename Type>
class polynomeBase : public containorSparse<Type>
{	
private:
	
protected:
	/**
	 * Méthode pour afficher une fonction de base
	 *
	 * Le problème est de pouvoir afficher tout les polynomes
	 * sans avoir à réécrire de méthode write ou d'opérateur <<.
	 * Pour cela on demande à toutes les classes enfants de réimplémenter
	 * cette méthode comme cela chaque classe affichera correctement
	 * la base de polynômes qui lui convient.
	 *
	 * Par exemple dans polynome le getBaseString affichera quelque
	 * chose comme X^d . Y^dy . Z^dz alors que pour un polynome centré
	 * on pourra afcher (X-CenterX)^dx . (Y-CenterY)^dy . (Z-CenterZ)^dz
	 *
	 * @param degree le degréé de la fonction de base du polynome dont 
	 * on souhaite connaitre la représentation.
	 */  
	virtual const string getBaseString(const int & degree) const;
  
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
	void add(const polynomeBase<Type> & p , const polynomeBase<Type> & q);

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
	void add(const polynomeBase<Type> & q);

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
	void divide(const polynomeBase<Type> & poly , const Type & scalar);

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
	void multiply(const polynomeBase<Type> & p , const polynomeBase<Type> & q);
	
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
	void multiply(const polynomeBase<Type> & q);

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
	void multiply(const polynomeBase<Type> & poly , const Type & scalar);

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
	 * Méthode de base pour faire la soustraction de deux polynomes.
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
	void soustract(const polynomeBase<Type> & p , const polynomeBase<Type> & q);

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
	void soustract(const polynomeBase<Type> & q);

public:
  
  /**
   * Constructeur par défaut.
   */
  polynomeBase(void);
  
  /**
   * Construcuteur de copie.
   *
	 * @param poly le polynome avec lequel le nouvel objet est initialisé.
	 */
	polynomeBase(const polynomeBase<Type> & poly);

	/**
   * Destucteur.
   */
  virtual ~polynomeBase(void);
    
  /**
   * Methode pour connaitre le degré du polynome.
	 *
	 * @return le degré du polynome.
   */
  int getDegree(void) const;
  
  /**
   * Mehode pour l'affchage du polynome.
	 *
	 * @param os le flux dans lequel on souahite écrire le polynome.
   */
  void write(ostream & os) const;
};


/**
 * Opérateur pour écrire un polynome dans un flux.
 * 
 * @param out le flux dans lequel on souhaite écrire le polynome.
 *
 * @param poly le polynome que l'on souahite écrire.
 */ 
template<typename Type>
ostream & operator<<(ostream & out , const polynomeBase<Type> & poly);


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Le constructeur par défaut.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
polynomeBase<Type>::polynomeBase(void)
: containorSparse<Type>()
{
 this->clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Le constructeur copie.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
polynomeBase<Type>::polynomeBase(const polynomeBase<Type> & poly)
: containorSparse<Type>()
{
  containorSparse<Type>::copy(poly);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Le destructeur.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
polynomeBase<Type>::~polynomeBase(void)
{
  this->clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire l'addition de deux polynomes.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase<Type>::add(const polynomeBase<Type> & p , const polynomeBase<Type> & q)
{
	containorSparse<Type>::copy(p);
	add(q);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire l'addition unaire de deux polynomes.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase<Type>::add(const polynomeBase<Type> & q)
{
	int degree;
	Type value;

	for(degree = q.begin() ; degree != q.end() ; degree = q.next(degree)) {
		value = q[degree] + this->getData(degree);
		setData(degree , value);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la division d'un polynome et d'un scalaire.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase<Type>::divide(const polynomeBase<Type> & poly , const Type & scalar)
{
	assert(scalar != 0);

	containorSparse<Type>::copy(poly);
	divide(scalar);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la division unaire d'un polynome et d'un scalaire.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase<Type>::divide(const Type & scalar)
{
	assert(scalar != 0);

	int i;
	Type value;
	
	for(i=(this->begin()) ; i != (this->end()) ; i= (this->next(i))) {
		value = this->getData(i) / scalar;
		setData(i,value);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour conaitre le degré du polynome.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
int polynomeBase<Type>::getDegree(void) const
{
	// !! Attention !! lorsque le polynome est vide : le polynome est nul.
	// Donc son degré est 0. Si on ne catch pas l'exception ici le rbegin
	// est un rend et cela renvoie -1. Ce qui est une abération.
	if(this->empty()) {
		return 0;
	}

	return this->rbegin();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour conaitre le degré du polynome.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
const string polynomeBase<Type>::getBaseString(const int & degree) const
{
	return monome2string("x",0,degree);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Methode pour faire la multiplication de deux polynomes
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase<Type>::multiply(const polynomeBase<Type> & p , const polynomeBase<Type> & q)
{
	int i,j;
	Type value;

	this->clear();

	for(i=q.begin() ; i != q.end() ; i=q.next(i)) {
		for(j=p.begin() ; j != p.end() ; j=p.next(j)) {
			value = this->getData(i+j) + q[i] * p[j];
			setData(i+j , value);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Methode pour faire la multiplication unaire de deux polynomes
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase<Type>::multiply(const polynomeBase<Type> & q)
{
	polynomeBase tmp(*this);
	multiply(tmp,q);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la multiplication unaire de d'un polynome et d'un scalaire.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase<Type>::multiply(const Type & scalar)
{
	int degree;
	Type value;

	if(scalar == 0) {
		this->clear();
		return;
	}
	
	for(degree = (this->begin()) ; degree != (this->end()) ; degree =(this->next(degree))) {
		value = scalar *  (this->getData(degree));
		setData(degree , value);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la multiplication de d'un polynome et d'un scalaire.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
 void polynomeBase<Type>::multiply(const polynomeBase<Type> & p , const Type & scalar)
{
	containorSparse<Type>::copy(p);
	multiply(scalar);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la soustraction de deux polynomes.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase<Type>::soustract(const polynomeBase<Type> & p , const polynomeBase<Type> & q)
{
	containorSparse<Type>::copy(p);
	soustract(q);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Méthode pour faire la soustraction unaire de deux polynomes.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase<Type>::soustract(const polynomeBase<Type> & q)
{
	int degree;
	Type value;

	for(degree = q.begin() ; degree != q.end() ; degree = q.next(degree)) {
		value = this->getData(degree) - q[degree];
		setData(degree , value);
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method WRITE for the polynome.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void polynomeBase<Type>::write(ostream & os) const
{
	bool has_previous = false;
	Type coefficient;
	int i;

	if(this->empty()) {
		os << "0.";
		return;
	}

	for( i=(this->begin()) ; i != (this->end()) ; i =(this->next(i))) {

		coefficient = this->getData(i);

		if(has_previous) {
			os << (coefficient > 0 ? " +" :" ");
		}

		has_previous = false;

		if(i==0) {
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
ostream & operator << (ostream & out , const polynomeBase<Type> & poly)
{
	poly.write(out);
	return out;
}

#endif
