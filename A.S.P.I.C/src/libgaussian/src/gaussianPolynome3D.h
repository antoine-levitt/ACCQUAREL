/* 
* The gaussian library of the A.S.P.I.C. 
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
#ifndef _GAUSSIAN_POLYNOME_3D_
#define _GAUSSIAN_POLYNOME_3D_

#include "gaussian3D.h"
#include <iomanip>
#include <iostream>
#include <polynome3D.h>
using namespace std;

class gaussianPolynome3D : public gaussian3D  , protected polynome3D 
{
private:

protected:
	/**
	 * Méthode pour faire la copie.
	 */
	void copy(const gaussianPolynome3D & gp);
	
	/**
	 * Méthode pour afficher la tete des monomes de base du polynome.
	 * et utiliser la méthode write de la clasee polynome correctement.
	 *
	 * @param i le degré du monome de base à afficher.
	 */
	const virtual string getBaseString(const ipoint<3> & i) const;
	
public:
	/**
	 * Le constructeur par défaut.
	 */
	gaussianPolynome3D(void);

	/**
	 * Le constructeur de copie.
	 */
	gaussianPolynome3D(const gaussianPolynome3D & gp);


	/**
	 * Constructeur avec une gaussienne.
	 */
	gaussianPolynome3D(const gaussian3D & g);

	/**
	 * Le destructeur.
	 */
	virtual ~gaussianPolynome3D(void);
		
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
	gaussianPolynome3D centerMonomeMultiply(const double & root , const int & dimension)const;

	/**	
	 * Methode qui multiplie un polynome par (variable-root)^d
	 * Autrement dit cela revient à ajouter root comme zéro du polynome.
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @param dimension la dimension pour laquelle on veut ajouter la racine.
	 *
	 * @return le polynome produit de (variable-root) et de l'objet.
	 *
	 * @warning la dimension ne peut prendre que les valeur 0,1 ou 2.
	 *
	 * @warning il faut que le degré soit positif ou nul.
	 */
	gaussianPolynome3D centerMonomeMultiply(const double & root , const int & dimension , const int & degree)const;

	/**
	 * Methode qui multiplie un polynome par (x-root)^degree
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @return le polynome produit de (x-root)^degree et de l'objet.
	 *
	 * @warning il faut que le degré soit positif ou nul.
	 */
	gaussianPolynome3D centerMonomeMultiplyX(const double & root , const int & degree =1)const;

	/**
	 * Methode qui multiplie un polynome par (y-root)^degree
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @return le polynome produit de (x-root)^degree et de l'objet.
	 *
	 * @warning il faut que le degré soit positif ou nul.
	 */
	gaussianPolynome3D centerMonomeMultiplyY(const double & root , const int & degree =1)const;

	/**
	 * Methode qui multiplie un polynome par (z-root)^degree
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @return le polynome produit de (z-root)^degree et de l'objet.
	 *
	 * @warning il faut que le degré soit positif ou nul.
	 */
	gaussianPolynome3D centerMonomeMultiplyZ(const double & root , const int & degree =1)const;

	/**	
	 * Methode qui multiplie un polynome par (variable-root)^d
	 * Autrement dit cela revient à ajouter root comme zéro du polynome.
	 *
	 * @param root la valeur du zéro à ajouter au polynome.
	 *
	 * @param dimension la dimension pour laquelle on veut ajouter la racine.
	 *
	 * @return le polynome produit de (variable-root) et de l'objet.
	 *
	 * @warning la dimension ne peut prendre que les valeur 0,1 ou 2.
	 *
	 * @warning il faut que le degré soit positif ou nul.
	 */
	gaussianPolynome3D centerMonomeMultiply(const dpoint<3> & root , const ipoint<3> & degree=ipoint<3>(1))const;

	/**
	 * Méthode pour changer le centre de la gaussienne.
	 *
	 * En fait la méthode setCenter permet de fixer la valeur du centre
	 * tout betement. Cette méthode à pour objectif de changer la valeur
	 * du centre et de conserver le polynome comme il était auparavant.
	 */
	void changeCenter(const dpoint<3> & center);

	/**
	 * Méthode pour mettre le polynome à 0.
	 */
	void clearPolynome(void);
	
		
	/**
	 * La dérivation.
	 *
	 * @param dim la direction x,y ou z dans laquelle on souhaite dériver.
	 *
	 * @warning le parametre dim ne peut prendre que les valeurs 0, 1 , ou 2.
	 */
	gaussianPolynome3D derivate(int dim) const;

	/**
	 * La dérivation suivant la variable z.
	 */
	gaussianPolynome3D derivateZ(void) const;

	/**
	 * La dérivation suivant la variable y.
	 */
	gaussianPolynome3D derivateY(void) const;

	/**
	 * La dérivation suivant la variable x.
	 */
	gaussianPolynome3D derivateX(void) const;

	/**
	 * Méthode générale de dérivation.
	 *
	 * @param degre un triplet d'entiers qui contient
	 * le nombre de dérivation à effectuer dans chaque direction.
	 */
	gaussianPolynome3D derivate(ipoint<3>  degree) const;

	/**
		* Method to eval the value of the gaussian polynome 
	 * in the point x.
	 */
	double eval(const dpoint<3> & position) const;
	
	/**
	 * Méthode pour accéder au polynome.
	 *
	 * @return le polynome qui se trouve devant la gaussienne.
	 */
	const polynome3D &  getPolynomeCoefficients(void) const;
	
	/**
	 * Méthode pour connaitre le coefficient d'un monome qui compose le polynome.
	 *
	 * @param degree le degré du monome dont on souhaite conniatre le coefficient.
	 *
	 * @return la valeur du coefficient du monome de degré degree.
	 *
	 * @warning il faut que toutes les composantes du triplet d'entier soient
	 * positives ou nulles.
	 */
	double getPolynomeCoefficient(const ipoint<3> & degree) const;

	/**
	 * Méthode pour connaitre le coefficient d'un monome qui compose le polynome.
	 *
	 * @param degreeX le degré en x du monome dont on souhaite conniatre le coefficient.
	 *
	 * @param degreeY le degré en y du monome dont on souhaite conniatre le coefficient.
	 *
	 * @param degreeZ le degré en z du monome dont on souhaite conniatre le coefficient.
	 *
	 * @return la valeur du coefficient du monome de degré degree.
	 *
	 * @warning il faut que le paramètre degreeX soit positif ou nul. 
	 *
	 * @warning il faut que le paramètre degreeZ soit positif ou nul. 
	 *
	 * @warning il faut que le paramètre degreeY soit positif ou nul. 
	 */
	double getPolynomeCoefficient(const int & degreeX , const int & degreeY , const int & degreeZ) const;

	/**
	 * Méthode pour connaitre les plus haut degré du polynome.
	 *
	 * @return le plus haut degré du polynome.
	 *
	 * @see polynome3D::getUniformMax4Degree.
	 */
	ipoint<3> getPolynomeUniformMax4Degree(void) const;

	/**
	 * méthode pour calculer le recouvrement.
	 *
	 * @return la valeur de l'intégrale sur R du polynome gaussien.
	 */
	double integral(void) const;

	/**
	 * La multiplication de deux gaussiennes polynomes.
	 */
	void multiply(gaussianPolynome3D gp_a , gaussianPolynome3D gp_b);

	/**
	 * La multiplication unaire de deux gaussiennes polynomes.
	 */
	void multiply(const gaussianPolynome3D & gaussianPoly);

	/**
	 * La multiplication d'une gaussienne polynome et d'in scalaire.
	 */
	void multiply(const gaussianPolynome3D & gaussianPoly , const double & scalar);

	/**
	 * La multiplication uanire d'une gaussienne polynome et d'un scalaire.
	 */
	void multiply(const double & scalar);

	/**
   * Opérateur pour faire la multiplication de deux gaussiennes
	 * polynome.
	 */
	gaussianPolynome3D operator * (const gaussianPolynome3D & gaussianPoly) const;

		/**
   * Opérateur pour faire la multiplication de deux gaussiennes
	 * polynome.
	 */
	gaussianPolynome3D operator* (const double & scalar) const;

	/**
   * Opérateur pour faire la multiplication unaire de deux gaussiennes
	 * polynome.
	 */
	gaussianPolynome3D & operator*= (const gaussianPolynome3D & gaussianPoly);

	/**
   * Opérateur pour faire la multiplication unaire de deux gaussiennes
	 * polynome.
	 */
	gaussianPolynome3D & operator*= (const double & scalar);

	/**
	 * Opérateur d'affectation avec un polynome gaussien.
	 *
	 * @param gaussianPoly le polynome gaussien à recopier dans l'objet.
	 *
	 * @return une référence vers l'objet modifié.
	 */
	gaussianPolynome3D & operator= (const gaussianPolynome3D & gaussianPoly);

	/**
	 * Opérateur d'affectation avec une gaussienne.
	 *
	 * @param g la gaussienne à recopier dans l'objet.
	 *
	 * @return une référence vers l'objet modifié.
	 */
	gaussianPolynome3D & operator= (const gaussian3D & g);

	/**
	 * Méthode pour connaitre le marqueur de fin du polynome
	 */
	ipoint<3> polynomeEnd(void) const;

	/**
	 * Méthode pour connaitre le premier degré du polynome
	 */
	ipoint<3> polynomeBegin(void) const;

	/**
	 * Méthode pour connaitre le suivant.
	 */
	ipoint<3> polynomeNext(const ipoint<3> & degree) const;

	/**
	 * Méthode pour donner au polynome une valeur constante.
	 *
	 * @param value la valeur constante que l'on souahite attribuer au polynome.
	 */
	void setPolynomeCoefficients(const double & value);

	/**
	 * Méthode pour donner une valeur au polynome à partir de coefficients.
	 *
	 * @param polynomeCoefficients les coefficients que l'on souhaite recopier dans
	 * l'objet.
	 */
	void setPolynomeCoefficients(const polynome3D & polynomeCoefficients);

	/**
	 * Méthode pour modifier la valeur d'un coefficient du polynome.
	 *
	 * @param degree le degree pour lequel on veut modifier le coefficient.
	 *
	 * @param value la valeur à attribuer au coefficient du polynome.
	 *
	 * @warning les trois composantes du triplet d'entier degree doivent etre
	 * positive ou nulles.
	 */
	void setPolynomeCoefficient(const ipoint<3> & degree , const double & value);

	/**
	 * Méthode pour fixer les degrés du polynome.
	 *
	 * @param degreeX le degré en x du monome dont on souhaite modifier le coefficient.
	 *
	 * @param degreeY le degré en y du monome dont on souhaite modifier le coefficient.
	 *
	 * @param degreeZ le degré en z du monome dont on souhaite modifier le coefficient.
	 *
	 * @param value la valuer que l'on souhaite attribuer au coefficient.
	 *
	 * @warning le parametre degreeX doit avoir une valeur positive ou nulle.
	 *
	 * @warning le parametre degreeY doit avoir une valeur positive ou nulle.
	 *
	 * @warning le parametre degreeZ doit avoir une valeur positive ou nulle.
	 */
	void setPolynomeCoefficient(const int & degreeX , const int & degreeY , const int & degreeZ , const double & value);
	
	/**
	 * Méthode pour écrire dans un flux.
	 *
	 * @param outStream le flux dans lequel on souhaite écrire l'objet.
	 */
	void write(ostream & outStream) const;
};

/**
* Opérateur pour écrire dans un flux.
*
* @param outStream le flux dans lequel on souhaite écrire l'objet.
*
* @param gaussianPoly le polynome gaussien que l'on souahite écrire dans le flux.
*/
extern ostream & operator<<( ostream & outStream , const gaussianPolynome3D & gaussianPoly);

#endif
