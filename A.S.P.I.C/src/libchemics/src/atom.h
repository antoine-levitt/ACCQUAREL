/* 
* The gaussian library of A.S.P.I.C. 
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
#ifndef _ATOM_
#define _ATOM_


#include "chemicalElement.h"
#include <dpoint.h>
#include <iostream>
#include <string>
using namespace std;

/**
* La classe "atom" est une classe pour la description et la manipulation des atomes.
 *
 * Au coeur de l'atome se trouve l'élément chimique qui lme compose. La classe atome hérite donc 
 * naturellement de cette classe. Cepedant l'atome, contrairement à l'élément chimique est 
 * distinguable d'un autre atome composé du même élément chimique par sa position et son nombre 
 * de neutrons.
 *
 * Un atome est donc un élément chimique pour lequel on a spécifie :
 * - Le nombre de neutrons
 * - La position
 * - L'unité dans laquelle on souhaite exprimer les distances.
 */
class atom : public chemicalElement
{
public:

	/**
	 * Enumération pour les unités de la position.
	 * En gros on peut utiliser des angströms ou des unités
	 * atomiques.
	 *
	 * - ANGSTROEM les distances sont exprimées en angströms (10E-10 m).
	 * - ATOMIC_UNIT les distances sont exprimées en unités atomiques et (1 a.u = rayon de l'atome de Bohre).
	 */
	typedef enum _distanceUnit {ANGSTROEM , ATOMIC_UNIT} distanceUnit;

private:
	/**
	 * L'unité pour les distances.
	 */
	distanceUnit DistanceUnit;

	/**
	 * Le nombre de neutrons de l'atome.
	 */
	int NbrOfNeutrons;

	/**
	 * La position de l'atome.
	 */
	dpoint<3> Position;
	
protected:

	/**
	 * Méthode pour copier un objet atome.
	 */
	void copy(const atom & a);

public :

	/**
	 * Constructeur par défaut.
	 */
	atom(void);
	
	/**
	 * Constructeur de copie.
	 *
	 * @param a l'atome a recopier dans l'objet nouvellement créé.
	 */
	 atom(const atom & a);

	 /**
	  * Destructeur de la classe.
		*/
	virtual ~atom(void);
	
	/**
	 * La valeur pour le rapport Unité atomique et Angström.
	 */
	static const double AtomicUnit2Angstroem;

	/**
	 * Méthode qui converti une distance exprimée en unité atomique
	 * en une distance exprimée en angström.
	 *
	 * @param la distance en angström.
	 *
	 * @return la meme distance en unité atomique.
	 */
	static double angstroem2AtomicUnitConverter(const double & distance);

	/**
	 * Méthode qui converti un triplet de distances (une position) exprimé en angström 
	 * en une position exprimée en unité atomique.
	 *
	 * @param la position en angström.
	 *
	 * @return la meme position en unité atomique.
	 */
	static dpoint<3> angstroem2AtomicUnitConverter(const dpoint<3> & position);

	/**
	 * Méthode qui converti une distance exprimée en unité atomique
	 * en une distance en angström.
	 *
	 * @param distance la distance en unité atomique.
	 *
	 * @return la meme distance en angström.
	 */
	static double atomicUnit2AngstroemConverter(const double & distance);

	/**
	 * Méthode qui converti un triplet de distances exprimées en unité atomique
	 * en une distance en angström.
	 *
	 * @param postion une position exprimée en unité atomique.
	 *
	 * @return la position en angström.
	 */
	static dpoint<3> atomicUnit2AngstroemConverter(const dpoint<3> & distance);

	 /**
	 * Méthode pour changer l'unité pour la position de l'atome.
	 *
	 * Cette méthode est "constante" car elle ne change pas l'objet
	 * à proprement parler même si elle change les valeurs de la 
	 * position de l'objet.
	 *
	 * @param unit l'unité dans laquelle on souahite exprimer la position.
	 */
	 void changeDistanceUnit(const distanceUnit & unit) const;

	/**
	 * Méthode pour remmettre un atome à "zéro".
	 */
	void clear(void);


	/**
	 * Méthode pour convertir une unité de distance en chaine de charctères.
	 *
	 * @param unit l'unité de distance que l'on souahite représenter par une chaine de charactères.
	 */
	static string distanceUnit2String(const distanceUnit & unit);

	/**
	 * Méthode GET pour l'unité de la position.
	 *
	 * @return l'unité de la position.
	 */
	const distanceUnit & getDistanceUnit(void) const;

	/**
	 * Méthode GET pour le nombre de neutrons.
	 *
	 * @return le nombre de neutrons de l'objet.
	 */
	const int & getNbrOfNeutrons(void) const;
	
	/**
	 * Méthode GET pour la position.
	 *
	 * @return la position de l'atome.
	 */
	const dpoint<3> & getPosition(void) const;

	/**
	 * Method GET for the atom position.
	 *
	 * @param unit the unit in which the position of the atom should be 
	 * expressed.
	 * 
	 * @return the position of the atom expressed in the unit unit.
	 */
	const dpoint<3> getPosition(const atom::distanceUnit & unit) const;
	
	/**
	 * Operateur d'affectation.
	 *
	 * @param a la valeur à affecter à l'objet.
	 */
	atom & operator=(const atom & a);
	
	/**
	 * MŽthode pour recopier la valeur d'un atome dans l'objet.
	 *
	 * @param a l'atome à recopier dans l'objet.
	 */
	void setAtom(const atom & a);
	 
	/**
	 * Méthode SET pour l'unité de la position.
	 *
	 * Cette méthode permet uniquement de spécifier l'unité. Pour les conversion
	 * il est fortement conseillé de se référeer à la méthode changePositionUnit.
	 * 
	 * @param unit l'unité pour la position.
	 */
	 void setDistanceUnit(const distanceUnit & unit);

	/**
	 * MŽthode pour modifier le nombre de neutrons.
	 * 
	 * @param nbrOfNeutrons le nombre de neutrons de l'objet.
	 *
	 * @warning le nombre de neutrons doit être positif ou nul.
	 */
	void setNbrOfNeutrons(const int & nbrOfNEutrons);
	
	/**
	 * Méthode pour modifier la position de l'atome.
	 *
	 * @param position la position de l'atome.
	 *
	 * @param unit l'unité dans laquelle on exprime la position.
	 */
	void setPosition(const dpoint<3> & position , const distanceUnit & unit=atom::ATOMIC_UNIT);

	/**
	 * Méthode pour modifier la position de l'atome.
	 *
	 * @param x la première coordonée cartésienne de la position.
	 *
	 * @param y la seconde coordonée cartésienne de la position.
	 *
	 * @param z la troisème coordonée cartésienne de la position.
	 *
	 * @param unit l'unité dans laquelle on exprime la poistion.
	 */
	 void setPosition(const double & x ,const double & y ,const double & z , const distanceUnit & unit=atom::ATOMIC_UNIT);
	
	/**
	 * Méthode pour convertir une chaîne de carctères en unité de distance.
	 *
	 * @param distanceUnitString une chaine de charactères sensée représenter une unité de distance.
	 * Les chaine connues sont : "au" , "AtomicUnit" , pour les untiés atomiques et 
	 * 
	 * @TODO.
	 *
	 * @return l'unité, lorsque celle ci à été trouvée.
	 */
	 static atom::distanceUnit string2DistanceUnit(const string & distanceUnitString);

	 /**
	 * Méthode pour écrire un atome dans un flux.
	 *
	 * @param outStream le flux dans lequel on souhaite écrire l'atome.
	 *
	 * @param linePrefix le préfixe que l'on souahite mettre avant chaque début de ligne.
	 */
	void write(ostream & outStream , const string & linePrefix = "") const; 


	/**
	 * Méthode pour écrire un atome dans un flux.
	 *
	 * @param outStream le flux dans lequel on souhaite écrire l'atome.
	 *
	 * @param linePrefix le préfixe que l'on souahite mettre avant chaque début de ligne.
	 */
	void writeXML(ostream & outStream , const string & linePrefix = "") const; 

};

/**
 * Operateur externe pour écrire un objet atome dans un flux.
 *
 * @param outStream le flux dans lequel on souhaite écrire l'atome.
 *
 * @param a l'atome que l'on souhaite décrire.
 *
 */
extern ostream & operator<<(ostream & outStream , const atom & a);


/**
 * Operateur externe pour écrire un objet atome dans un flux.
 *
 * @param outStream le flux dans lequel on souhaite écrire l'atome.
 *
 * @param a l'atome que l'on souhaite décrire.
 *
 */
extern ostream & operator<<(ostream & outStream , const atom::distanceUnit & unit);


#endif

