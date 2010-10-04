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
#ifndef POLYNOME3D
#define POLYNOME3D

#include "polynome.h"
#include "polynomeBase3D.h"

class polynome3D : public polynomeBase3D<double>
{
	/**
	 * Methode pour afficher les fonction de base.
	 *
	 * Comme nous travaillons ici avec les polynômes usuel la 
	 * fonction de base pour un degré d, s'écriit "x^d".
	 *
	 * @param i le degré de la fonction de base qui l'on cherche à convertir en
	 * chaine de charactères.
	 *
	 * @return la chaine de charactère qui repésente la fonction de base, ie "1" lorsque
	 * le degré est nul et "x^d" lorsque le degré demandé est supérieur ou égal à 1.
	 */
	virtual const string  getBaseString(const ipoint<3> & i) const;

public:

	/**
	 * Constructeur par defaut.
	 */
	polynome3D(void);

	/**
	 * Constructeur de copie.
	 *
	 * @param poly le polynome avec lequel on souhaite initialiser l'objet.
	 */
	polynome3D(const polynome3D & poly);

	/**
	 * Construcuteur de polynome constant.
	 *
	 * Ce constructeur permet de créer un polynome constant dont la
	 * valeur est passée en argument.
	 */
	polynome3D(const double & value);

	/**
	 * Destructeur.
	 */
	virtual ~polynome3D(void);

	/**
	 * Methode qui multiplie un polynome par (variable-root)
	 * Autrement dit cela revient à ajouter root comme zéro du polynome.
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @param dimension la dimension pour laquelle on veut ajouter la racine.
	 *
	 * @return le polynome produit de (variable-root) et de l'objet.
	 *
	 * @warning la dimension ne peut prendre que les valeur 0,1 ou 2.
	 */
	polynome3D centerMonomeMultiply(const double & root , const int & dimension)const;

	/**	
	 * Methode qui multiplie un polynome par (variable-root)^d
	 * Autrement dit cela revient à ajouter root comme zéro du polynome.
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @param dimension la dimension pour laquelle on veut ajouter la racine.
	 *
	 * @param degree l'ordre de multiplicité de la racine à ajouter.
	 *
	 * @return le polynome produit de (variable-root) et de l'objet.
	 *
	 * @warning la dimension ne peut prendre que les valeur 0,1 ou 2.
	 *
	 * @warning il faut que le degré soit positif ou nul.
	 */
	polynome3D centerMonomeMultiply(const double & root , const int & dimension , const int & degree)const;

	/**
	 * Methode qui multiplie un polynome par (x-root)^degree
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @param degree l'ordre de multiplicité de la racine à ajouter.
	 *
	 * @return le polynome produit de (x-root)^degree et de l'objet.
	 *
	 * @warning il faut que le degré soit positif ou nul.
	 */
	polynome3D centerMonomeMultiplyX(const double & root , const int & degree =1)const;

	/**
	 * Methode qui multiplie un polynome par (y-root)^degree
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @param degree l'ordre de multiplicité de la racine à ajouter.
	 *
	 * @return le polynome produit de (x-root)^degree et de l'objet.
	 *
	 * @warning il faut que le degré soit positif ou nul.
	 */
	polynome3D centerMonomeMultiplyY(const double & root , const int & degree =1)const;

	/**
	 * Methode qui multiplie un polynome par (z-root)^degree
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @param degree l'ordre de multiplicité de la racine à ajouter.
	 *
	 * @return le polynome produit de (z-root)^degree et de l'objet.
	 *
	 * @warning il faut que le degré soit positif ou nul.
	 */
	polynome3D centerMonomeMultiplyZ(const double & root , const int & degree =1)const;

	/**	
	 * Methode qui multiplie un polynome par (variable-root)^d
	 * Autrement dit cela revient à ajouter root comme zéro du polynome.
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @param dimension la dimension pour laquelle on veut ajouter la racine.
	 *
	 * @param degree l'ordre de multiplicité de la racine à ajouter.
	 *
	 * @return le polynome produit de (variable-root) et de l'objet.
	 *
	 * @warning la dimension ne peut prendre que les valeur 0,1 ou 2.
	 *
	 * @warning il faut que le degré soit positif ou nul.
	 */
	polynome3D centerMonomeMultiply(const dpoint<3> & root , const ipoint<3> & degree=ipoint<3>(1))const;

	/**
	 * Methode de base pour la dérivation.
	 *
	 * Cette méthode permet de dériver un polynome par rapport à la variable x , y ou z
	 * suivant la valeur du parametre dimension.
	 *
	 * @param dimension dimension par rapport à laquelle on souhaite dériver le polynome, 
	 * lorsque dimension prend la valeur 0 on dérive par rapport à , 1 par rapport à y et 2
	 * par rapport à z.
	 *
	 * @return le polynome dérivé.
	 *
	 * @warning la vaelrude la dimension ne peut etre que 0, 1 ou 2.
	 */
	polynome3D derivate(const int & dimension) const;

	/**
	 * Methode qui dérive un polynome par rapport à la variable x.
	 *
	 * @param derivationDegree le degré de dérivation du polynome, sachant que
	 * la première composante du triplet d'entier contient la valeur du nombre de dérivation
	 * par rapport à la première variable, la seconde la valeur du nombre de dérivation par 
	 * rapport à la seconde variable, et la derniere le nombre de dérivation par rapport
	 * à la dernière variable.
	 *
	 * @return le polynome issu de la dérivation.
	 */
	polynome3D derivate(const ipoint<3> & derivationDegree) const;

	/**
	 * Methode qui dérive un polynome par rapport à la variable x.
	 *
	 * @return le polynome dérivé par rapport à la première variable.
	 */
	polynome3D derivateX(void) const;
	
	/**
	 * Methode qui dérive un polynome par rapport à la variable y.
	 *
	 * @return le polynome dérivé par rapport à la seconde variable.
	 */
	polynome3D derivateY(void) const;

	/**
	 * Methode qui dérive un polynome par rapport à la variable z.
	 *
	 * @return le polynome dérivé par rapport à la troisième variable.
	 */
	polynome3D derivateZ(void) const;
	
	
	/**
	 * Method to eval the value of the polynome 
	 * in the point x.
	 */
	double eval(const dpoint<3> & position) const;
	
	

	/**
	 * Methode qui multiplie un polynome par varaible^degree
	 *
	 * @param degree le degree du monome avec lequel souhaite multiplier le polynome.
	 *
	 * @param dimension la variable avec laquelle  on souahite opérer.
	 *
	 * @return le polynome qui contient le produit de x^d et de l'objet.
	 */
	polynome3D monomeMultiply(const int & dimension , const int & degree = 1) const;

	/**
	 * Methode qui multiplie un polynome par x^degreeX
	 *
	 * @param degree le degree du monome avec lequel souhaite multiplier le polynome.
	 *
	 * @return le polynome qui contient le produit de x^d et de l'objet.
	 */
	polynome3D monomeMultiplyX(const int & degree = 1) const;

	/**
	 * Methode qui multiplie un polynome par y^degree
	 *
	 * @param degree le degree du monome avec lequel souhaite multiplier le polynome.
	 *
	 * @return le polynome qui contient le produit de y^degree et de l'objet.
	 */
	polynome3D monomeMultiplyY(const int & degree = 1) const;
	
	/**
	 * Methode qui multiplie un polynome par z^degree
	 *
	 * @param degree le degree du monome avec lequel souhaite multiplier le polynome.
	 *
	 * @return le polynome qui contient le produit de z^degree et de l'objet.
	 */
	polynome3D monomeMultiplyZ(const int & degree = 1) const;
	
	/**
	 * Methode qui multiplie un polynome par x^degreeX * y^degreeY * z^degreeZ.
	 *
	 * @param degree le degree du monome avec lequel souhaite multiplier le polynome.
	 *
	 * @return le polynome qui contient le produit de x^d et de l'objet.
	 */
	polynome3D monomeMultiply(const ipoint<3> & degree = ipoint<3>(1)) const;

	/**
	 * Operateur pour l'assignement.
	 * 
	 * @param p le polynome que l'on souhaite copier dans
	 * l'objet appelant.
	 *
	 * @return une référence vers le nouveau polynome.
	 */
	polynome3D & operator= (const polynome3D & p);

	/**
	 * Operateur pour faire l'addition de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite ajouter à l'objet.
	 *
	 * @return un polynome qui contient la somme demandée.
	 */
	polynome3D operator+ (const polynome3D & p) const;

	/**
	 * Opérateur pour faire une addition unaire de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite ajouter à l'objet.
	 *
	 * @return une référence vers le polynome qui contient la somme demandée.
	 */
	polynome3D & operator+= (const polynome3D & p);	
	

	/**
	 * Operateur pour faire la soustraction de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite soustraire ajouter à l'objet.
	 *
	 * @return un polynome qui contient la différence demandée.
	 */
	polynome3D operator- ( const polynome3D & p) const;

	/**
	 * Operateur pour faire la soustraction unaire de deux polynomes.
	 *
	 * @param p le polynome que l'on souhaite soustraire à l'objet.
	 *
	 * @return une référence vers le polynome qui contient la différence demandée.
	 */
	polynome3D & operator-= (const polynome3D & p);

	/**
	 * Operateur pour faire la multiplication de deux polynomes.
	 * 
	 * @param p le polynome avec lequel on souhaite multiplier l'objet.
	 *
	 * @return un polynome qui contient le produit demandé.
	 */
	polynome3D operator* (const polynome3D & p) const;

	/**
	 * Operateur pour faire la multiplication d'un polynome et d'un scalaire.
	 * 
	 * @param scalar le scalaire avec lequel on souhaite multiplier l'objet.
	 *
	 * @return un polynome qui contient le produit demandé.
	 */
	polynome3D operator* (const double & scalar) const;

	/**
	 * Operateur pour faire une multiplication unaire de deux polynomes.
	 * 
	 * @param p le polynome avec lequel on souhaite multiplier l'objet.
	 *
	 * @return une référence vers le polynome qui contient le produit demandé.
	 */
	polynome3D & operator*= (const polynome3D & p);

	/**
	 * Operateur pour faire une multiplication unaire d'un polynome et d'un scalaire.
	 * 
	 * @param scalar le scalaire avec lequel on souhaite multiplier l'objet.
	 *
	 * @return une référence vers le polynome qui contient le produit demandé.
	 */
	polynome3D & operator*= (const double & scalar);

	/**
	 * Operateur pour faire la disision d'un polynome et d'un scalaire.	
	 *
	 * @param scalar le scalaire par lequel on souhaite diviser l'objet.
	 *
	 * @return un polynome qui contient le quotient demandé.
	 */
	polynome3D operator/ (const double & scalar) const;

	/**
	 * Operateur pour faire une division unaire d'un polynome et d'un sclaire.
	 *
	 * @param scalar le scalaire par lequel on souhaite diviser l'objet.
	 *
	 * @return une référence vers le polynome qui contient le quotient demandé.
	 */
	polynome3D & operator/= (const double & scalar);

	/**
	 * Méthode pour changer la varialble linéaire.
	 *
	 * @param scalar le facteur du changement de variable.
	 *
	 * @return la forme du polynome après le changement de variable.
	 */
	polynome3D scale(const dpoint<3> & scalar);

	/**
	 * Méthode pour donner au polynome une valeur constante.
	 *
	 * @param value la valeur que l'on souhaite donner au polynome.
	 */
	void setPolynomeCoefficients(const double & value);

	/**
	 * Méthode pour recopier les coefficients d'un polynome.
	 *
	 * @param poly les coefficients à recopier dans l'objet.
	 */
	void setPolynomeCoefficients(const polynome3D & poly);

	/**
	 * Méthode pour construire un polynome tri dimensionel comme un produit de
	 * polynome mono dimensionnel.
	 *
	 * @param polyX le polynome en x.
	 * @param polyY le polynome en y.
	 * @param polyZ le polynome en z.
	 */
	void setPolynomeCoefficients(const polynome & polyX , const polynome & polyY , const polynome & polyZ);

};

////////////////////////////////////////////////////////////////////////////////
// Operateurs externe de multiplication pour les scalaires :
//
// 1- La multplication :
extern polynome3D operator*(const double & scalar , const polynome3D & p);


#endif

