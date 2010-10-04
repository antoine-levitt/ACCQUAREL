#ifndef _XML_BASIS_KEY_LIST_PARSER_
#define _XML_BASIS_KEY_LIST_PARSER_

#include "aspicConfiguration.h"
#include <containor.h>
#include <xmlDocumentParser.h>

/**
 * Classe pour retrouver une liste de toutes les clés de base disponible sur 
 * la plate forme.
 */
class xmlBasisKeyListParser : public xmlDocumentParser
{
private:
	
	/**
	 * L'objet xmlBasisKeyListParser qui va vraiment tout faire.
	 */
	static xmlBasisKeyListParser * BasisKeyListParser;

protected:

	/**
	 * Constructeur par défaut.
	 */
	inline xmlBasisKeyListParser(void)
	{
		;
	}

	/**
	 * Destructeur.
	 */
	inline virtual ~xmlBasisKeyListParser(void)
	{
		;
	}

	/**
	 * Méthode pour retrouver le schéma associé au fichier XML que nous allons parsé.
	 */
	inline string getSchemaURI(void) const
	{
		return aspicConfiguration::getBasisKeyListSchemaPath();
	}


	/**
	 * Méthode pour retrouver dans l'objet la liste des clés pour les bases gaussiennes.
	 *
	 */
	inline const containor<string> getBasisList(void) const
	{
		containor<string> basisList;
		DOMNodeList * nodeBasisList;
		int nbrOfBasis , i;

		nodeBasisList = getElementsByTagName("BasisKey");
		
		if(nodeBasisList == NULL || (nbrOfBasis=nodeBasisList->getLength()) == 0) {
			basisList.clear();
			return basisList;
		} 
			
		basisList.setSizes(nbrOfBasis);

		for(i=0 ; i < nbrOfBasis ; i++) {
			basisList[i] = getNodeStringValue(nodeBasisList->item(i),xmlParser::Remove_White_Space);
		}

		return basisList;
	}

public:

	/**
	 * Méthode pour se connecter.
	 */
	inline static void connect(void) {
	
		if(connected()) {
			return;
		}
	
		BasisKeyListParser = new xmlBasisKeyListParser;
		BasisKeyListParser->load(aspicConfiguration::getBasisKeyListPath());
	}

	/**
	 * Méthode pour savoir si l'on est déconnecté ou pas.
	 *
	 * @return true lorsque l'on est connecté, false sinon.
	 */
	inline static bool connected(void)
	{
		if(BasisKeyListParser == NULL) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * Méthode pour se déconnecter.
	 */
	inline static void disconnect(void)
	{
		if(connected() == false) {
			return;
		} else {
			delete BasisKeyListParser;
			BasisKeyListParser = NULL;
		}
	}


	/**
	 * Méthode pour récupérer la liste des bases gaussiennes présentes :
	 *
	 * @return la liste de toutes les clés de base qui peuvent etre utilisées
	 * sur la plate forme.
	 */
	static const containor<string> getBasisKeyList(void)
	{
		containor<string> basisKeyList;

		if(connected() == false) {
			connect();
		}	 
		
		basisKeyList = BasisKeyListParser->getBasisList();

		disconnect();
		return basisKeyList;
	}

};

#endif
