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
#ifndef _XML_CONTRACTION_PARSER_
#define _XML_CONTRACTION_PARSER_


#include "contractions.h"
#include <xmlParser.h>

class xmlContractionsParser : public xmlParser
{

private:
			
protected:

public:
	/**
	 * Constrcuteur par défaut.
	 */
	xmlContractionsParser(void);

	/**
	 * Constructeur.
	 */
	xmlContractionsParser(DOMElement * xmlContractionsElement);

	/**
	 * Destructeur.
	 */
	virtual ~xmlContractionsParser(void);
	
	/**
	 * Method GET for all contractions.
	 */
	contractions getContractions(void) const;

		/**
	 * Method GET for all contractions.
	 */
	contractions getContractions(const DOMElement * rootElement) const;

	/**
	 * Méthode pour retrouver le nom de la liste des contractions.
	 */
	static const string getContractionsListTagName(void);

	/**
	 * The name of the tag containing contraction pair data.
	 */
	static const string getContractionTagName(void);
	
	/**
		* The name of the tag containing the exponent data.
	 */
	static const string getExponentTagName(void);
	
	/**
		* The name of the tag containing the coefficient data.
	 */
	static const string getCoefficientTagName(void);
	
};

#endif
