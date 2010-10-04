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
#ifndef _MOLECULE_DOCUMENT_PARSER_
#define _MOLECULE_DOCUMENT_PARSER_

#include "molecule.h"
#include <xmlDocumentParser.h>
#include "xmlMoleculeParser.h"

class moleculeDocumentParser : public xmlDocumentParser , public xmlMoleculeParser
{
	private:

	protected:

	/**
	 * Method for the schema location.
	 */
	virtual string getSchemaURI(void) const;

	public:

	/**
	 * The constructor.
	 */
	moleculeDocumentParser(void);

	/**
	 * the destructor.
	 */
	~moleculeDocumentParser(void);

	/**
	 * Method GET for a molecule.
	 */
	molecule getMolecule4Document(const string & fileName);
};

#endif

