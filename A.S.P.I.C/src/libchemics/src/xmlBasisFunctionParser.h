/* 
 * The chemics library of A.S.P.I.C. 
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
#ifndef _GAUSSIAN_BASIS_FUNCTION_PARSER_
#define _GAUSSIAN_BASIS_FUNCTION_PARSER_

#include "atom.h"
#include "gaussianBasisFunction.h"
#include<dpoint.h>
#include <iostream>
#include <ipoint.h>
#include <string>
#include <xmlParser.h>
using namespace std;



/**
 * Classe pour lire la description d'une fonction de base dans 
 * un fichier XML.
 *
 * Un exemple de fichier est le suivant :
 *
 * <BasisFunction>
 *
 *		<Position>
 *			<x> 0 </x>
 *			<y> 0 </y>
 *			<z> 0 </z>
 *		</Position>
 *	
 *	<MonomeDegree>
 *		<x> 0 </x>
 *		<y> 0 </y>
 *		<z> 0 </z>
 *	</MonomeDegree>
 *
 *	<ContractionsList>
 *		<Contraction>
 *			<Coefficient>0.2149354488921539</Coefficient>
 *			<Exponent>18.73113696</Exponent>
 *		</Contraction>
 *		<Contraction>
 *			<Coefficient>0.3645712021918747</Coefficient>
 *			<Exponent>2.825394365</Exponent>
 *		</Contraction>
 *		<Contraction>
 *			<Coefficient>0.4150514278828989</Coefficient>
 *			<Exponent>0.6401216923</Exponent>
 *		</Contraction>
 *	</ContractionsList>
 *</BasisFunction>
 *
 * Cela permet de décrire :
 * - <Center> Le centre de la fonction de base (postion + unité)
 * - <MonomeDegree> le degré du polynome de la fonction de base.
 * - <ContractionList> la liste des contractions qui conposent la fonction de base.
 */
class xmlBasisFunctionParser : public virtual xmlParser
{
private:
protected:
public:

	/**
	 * Default constructor.
	 */
	xmlBasisFunctionParser(void);

	/**
	 * Constructor.
	 */
	xmlBasisFunctionParser(DOMElement * rootElement);
	
	/**
	 * Destructor.
	 */
	virtual ~xmlBasisFunctionParser(void);
	
	/**
	 * Method GET for a gaussian basis function.
	 */
	gaussianBasisFunction getBasisFunction(void) const; 

	/**
	 * Method GET for a gaussian basis function.
	 */
	gaussianBasisFunction getBasisFunction(const DOMElement * rootElement) const; 

	/**
	 * Méthode pour retrouver le nom du tag qui contient la fonction de base.
	 */
	const string getBasisFunctionTagName(void) const;

	/**
	 * Méthode pour retrouver le nom du tag qui contient le degré 
	 * du monome de la fonction de base.
	 */
 const string getMonomeDegreeTagName(void) const;

 /**
	 * Méthode pour retrouver le nom du tag qui contient le degré 
	 * du monome de la fonction de base.
	 */
 const string getCenterTagName(void) const;

};

#endif

