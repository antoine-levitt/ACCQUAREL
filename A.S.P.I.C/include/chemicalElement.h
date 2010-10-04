/* 
* The chemics library of A.S.P.I.C. 
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
#ifndef _CHEMICAL_ELEMENT_
#define _CHEMICAL_ELEMENT_

#include <assert.h>
#include <containor.h>
#include <iostream>
#include <map>
#include <string>
using namespace std;

/**
 *  Le problème des conteneurs , c'est qu'il veulent pouvoir afficher. Ici on 
 * leur fourni une méthode.
 */
extern ostream & operator<<(ostream & outStream , const pair<int,double> value);

/**
 * A chemical element is typically an object that describes structure of
 * a particualr element.
 */
class chemicalElement 
{
private:
	
	/**
	 * Un tableau qui contient les isotopes de l'élément chimique.
	 *
	 * En fait nous avons un tableau avec les paires (nombre de neutrons, probabilité)
	 */
	containor<pair<int,double> > Isotopes;
	
	/**
	 * The key of the chemical Element in the data base.
	 */
	string Key;

	/**
	 * the name of the chemical element. 
	 * For example Hydrogen.
	 */
	string Name;
	
	/**
	 * The number of protons of the chemical element.
	 * For example 1.
	 */
	int NbrOfProtons;
	
	/**
	 * the symbol of the chemical element.
	 * For Example H.
	 */
	string Symbol;

protected:

		/**
		 * Method that copys a chemical element.
		 *
		 * @param element the chemical element that will be copy in
		 * the calling oject.
		 */
		void copy(const chemicalElement & element);
	
public:

		/**
		 * Constructor.
		 * Constructs an empty chemical element.
		 */
		chemicalElement(void);
	
		/**
		 * Constructor.
		 * The copy constructor.
		 *
		 * @param element the element the calling object
		 * is initialize with.
		 */
		chemicalElement(const chemicalElement & element);
		
		/**
		 * Destructor.
		 */
		~chemicalElement(void);
		
		/**
     * Method that clears achemicalElement.
		 */
		void clear(void);

		/**
		 * Méthode pour savoir si l'objet manipulé est vide.
		 *
		 * @return true lorsque l'objet manipulé est vide et 
		 * false sinon.
		 */
		bool empty(void) const;

		/**
		 * Method GET for the key.
		 */
		const string & getChemicalElementKey(void) const;

		/**
		 * Méthode get pour connaitre le nombre d'isotopes connus.
		 */
		const pair<int , double> & getIsotope(const int & item) const;

		/**
		 * Méthode get pour connaitre le nombre d'isotopes connus.
		 */
		const containor<pair<int , double> > & getIsotopes(void) const;

		/**
		 * Method GET for the Name.
		 *
		 * @return the name of the chemical element.
		 */
		const string & getName(void) const;
		
		/**
		 * Méthode get pour connaitre le nombre d'isotopes connus.
		 */
		const int & getNbrOfIsotopes(void) const;

		/**
		 * Méthode get pour connaitre le nombre de neutrons d'un isotope.
		 */
		const int & getNbrOfNeutrons4Isotope(const int & item) const;

		/**
		 * Méthode get pour connaitre la probabilité d'existance d'un isotope.
		 */
		const double & getProbability4Isotope(const int & item) const;

		/**
		 * Method GET for the NbrOfProtons.
		 *
		 * @return the number of protons of the chemical element.
		 */
		const int & getNbrOfProtons(void) const;
		
		/**
		 * Method GET for the Symbol.
		 *
		 * @return the symbol of the chemical element.
		 */
		const string & getSymbol(void) const;
		
		/**
		 * Affectation operator.
		 * This enable the "=" formalism with chemicalElements.
		 *
		 * @param element the element to put in the calling object.
		 */
		chemicalElement & operator=(const chemicalElement & element);
		
		/**
	   * Method SET for a chemical Element.
	   */
		void setChemicalElement(const chemicalElement & element);
		
		/**
		 * Method SET for a chemical Element.
		 */
		void setChemicalElement(const string & chemicalElementKey);
		
		/**
		 * Method GET for the key.
		 */
		void setChemicalElementKey(const string & chemicalElementKey);

		/**
		 * Méthode pour fixer les isotopes à partir d'un tableau.
		 */
			void setIsotope(const int & item , const pair<int,double> & isotope);
	
		/**
		 * Méthode pour fixer les isotopes à partir d'un tableau.
		 */
		void setIsotope(const int & item , const int & nbrOfNeutrons , const double & probability);
	
		/**
		 * Méthode pour fixer les isotopes à partir d'un tableau.
		 */
		void setIsotopes(const containor<pair<int , double> > & isotopes);

		/**
		 * Method SET for the Name.
		 *
		 * @param name the name for the chemical element.
		 */
		void setName(const string & name);
		
		/**
		 * Méthode pour fixer le nombre d'isotopes de l'élément chimique.
		 */
		void setNbrOfIsotopes(const int & nbrOfIsotopes);

		/**
		 * Method SET for the NbrOfProtons.
		 *
		 * @param nbrOfPronts the number of protons of the chemical element.
		 *
		 * @warning the number of proton in argument must be non negative.
		 */
		void setNbrOfProtons(const int & nbrOfProtons);
		
		/**
		 * Method SET for the Symbol.
		 *
		 * @param symbol the symbol of the chemical element.
		 */
		void setSymbol(const string & symbol);

		/**
		 * Method to write a chemical element in a stream.
		 */
		void write(ostream & out , const string & prefixLine = "") const;
};

/**
 * extern operator << 
 */
ostream & operator<<(ostream & out, const chemicalElement & element);
	
#endif


