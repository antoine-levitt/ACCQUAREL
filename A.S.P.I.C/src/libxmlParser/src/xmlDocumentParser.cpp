/* 
* The xml parser library of A.S.P.I.C. 
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
#include "xmlDocumentParser.h"


///////////////////////////////////////////////////////////////
// Constructeur.
///////////////////////////////////////////////////////////////
xmlDocumentParser::xmlDocumentParser(void)
	: xmlParser() , 
	Parser(NULL) , 
	Document(NULL)
{
	try {
		XMLPlatformUtils::Initialize();
	} catch(const XMLException &toCatch) {
		cerr << "Error : in gaussianBasisFunctionParser::gaussianBasisFunctionParser(void)." << endl;
		cerr << "Error : an exception occurs during Xerces-c Initialization." << endl;
		char * message = XMLString::transcode(toCatch.getMessage());
		cerr << "Error :  exception message :" << message << endl;
   	XMLString::release(&message);
		exit(1);
	}
}

///////////////////////////////////////////////////////////////
// Destructeur.
///////////////////////////////////////////////////////////////
xmlDocumentParser::~xmlDocumentParser(void)
{
	close();

	// And call the termination method
  XMLPlatformUtils::Terminate();

}

///////////////////////////////////////////////////////////////
// Destructeur.
///////////////////////////////////////////////////////////////
void xmlDocumentParser::close(void)
{
	// On va quand même détruire le parser.
	if(Parser != NULL) {
		delete Parser;
		Parser = NULL;
		Document = NULL;
	}

	setRootElement(NULL);

}

///////////////////////////////////////////////////////////////
// Destructeur.
///////////////////////////////////////////////////////////////
const DOMDocument * xmlDocumentParser::getDocument(void) const
{
	return Document;
}

void xmlDocumentParser::load(const string & fileName)
{
	////////////////////////////////////////////////////////////
	// Création du parser Cerces.
	////////////////////////////////////////////////////////////
	if(Parser != NULL) {
		delete Parser;
	}

	Parser = new XercesDOMParser;

	/////////////////////////////////////////////////////////////
	// - Option : setDoNamespaces
	// true : Perform Namespace processing.
	//
	// false : Do not perform Namespace processing. 
	//
	// Default : false.
	// 
	// Note :  If the validation scheme is set to Val_Always 
	// or Val_Auto, then the document must contain a grammar 
	// that supports the use of namespaces.  
	/////////////////////////////////////////////////////////////
	Parser->setDoNamespaces(true);

	/////////////////////////////////////////////////////////////
	// - Option : setDoSchema
	//
	// true : Enable the parser's schema support.
	// 
	// false :  Disable the parser's schema support.
	//
	// Default : false.
	// 
	// Note :  If set to true, namespace processing must also be turned on.  
	/////////////////////////////////////////////////////////////
	Parser->setDoSchema(true);
	
	/////////////////////////////////////////////////////////////
	// - Option : setValidationSchemaFullChecking
	//
	// true :  Enable full schema constraint checking, including 
	// checking which may be time-consuming or memory intensive. 
	// Currently, particle unique attribution constraint checking 
	// and particle derivation restriction checking are controlled 
	// by this option.
	//
	// false : Disable full schema constraint checking.
	// 
	// Default : false
	//
	// Note :  This feature checks the Schema grammar itself for 
	// additional errors that are time-consuming or memory 
	// intensive. It does not affect the level of checking 
	// performed on document instances that use Schema grammars. 
	/////////////////////////////////////////////////////////////
	Parser->setValidationSchemaFullChecking(false);

	//////////////////////////////////////////////////////////////
	// - Option : setCreateEntityReferenceNodes
	// 
	// true :  Create EntityReference nodes in the DOM tree. 
	// The EntityReference nodes and their child nodes will be read-only. 
	//
	// false : Do not create EntityReference nodes in the DOM tree. 
	// No EntityReference nodes will be created, only the nodes 
	// corresponding to their fully expanded substitution text will 
	// be created.
	// 
	// Default : true.
	//
	// Note :  This feature only affects the appearance of 
	// EntityReference nodes in the DOM tree. The document will 
	// always contain the entity reference child nodes.  
	//////////////////////////////////////////////////////////////
	Parser->setCreateEntityReferenceNodes(true);

	///////////////////////////////////////////////////////////////
	// - Option : setIncludeIgnorableWhitespace(const bool) 
	//
	// true : Include text nodes that can be considered 
	// "ignorable whitespace" in the DOM tree.  
	//
	// false : Do not include ignorable whitespace in the DOM tree.
	//
	// Default : true
	//
	// Note :  The only way that the parser can determine if text 
	// is ignorable is by reading the associated grammar and having 
	// a content model for the document. When ignorable whitespace text 
	// nodes are included in the DOM tree, they will be flagged as 
	// ignorable; and the method DOMText::isIgnorableWhitespace() 
	// will return true for those text nodes.
	////////////////////////////////////////////////////////////////
	Parser->setIncludeIgnorableWhitespace(false); 
	
	/////////////////////////////////////////////////////////////////
	//  setExternalNoNamespaceSchemaLocation(const XMLCh* const) 
	//
	// Description : The XML Schema Recommendation explicitly 
	// states that the inclusion of schemaLocation/ 
	// noNamespaceSchemaLocation attributes in the instance document 
	// is only a hint; it does not mandate that these attributes must 
	// be used to locate schemas. This property allows the user to specify 
	// the no target namespace XML Schema Location externally. If 
	// specified, the instance document's noNamespaceSchemaLocation 
	// attribute will be effectively ignored. 
	//
	// Value : The syntax is the same as for the noNamespaceSchemaLocation 
	// attribute that may occur in an instance document: e.g."file_name.xsd"
	//////////////////////////////////////////////////////////////////
	string schemaFileName = getSchemaURI();
	
	if(!schemaFileName.empty()) {
		/////////////////////////////////////////////////////////////
	// Option : setValidationScheme.
	//
	// Val_Auto:  The parser will report validation errors only if a grammar is
	// specified. 
	//
	// Val_Always:  The parser will always report validation errors.
	//
	// Val_Never: Do not report validation errors.  
	//
	// default:  Val_Auto 
	//
	// Note :  If set to Val_Always, the document must specify a grammar. 
	// If this feature is set to Val_Never and document specifies a grammar,
	// that grammar might be parsed but no validation of the document 
	// contents will be performed. 
	/////////////////////////////////////////////////////////////
	Parser->setValidationScheme(XercesDOMParser::Val_Always);

		setExternalNoNamespaceSchemaLocation(schemaFileName);
	}

	///////////////////////////////////////////////////////////////
	// - Rapport d'erreur.
	///////////////////////////////////////////////////////////////
	xmlErrorReporter * errorReporter = new xmlErrorReporter;
	Parser->setErrorHandler(errorReporter);


	/////////////////////////////////////////////////////////////////
	// Parsage du fichier.
	/////////////////////////////////////////////////////////////////
	bool exceptionOccurs = false;
	try {
		Parser->parse(fileName.c_str());
	} catch (const OutOfMemoryException & ) {
		cerr << "Error : in gaussianBasisFunction::load(void)" << endl; 
		cerr << "Error : out of memory exception" << endl;
	 	exceptionOccurs = true;
	
	} catch(const XMLException & toCatch) {
		
		/////////////////////////////////////////////////////////////////////////////
		// Gestion des exceptions XML.
		/////////////////////////////////////////////////////////////////////////////
		cerr << "Error : in gaussianBasisFunction::load(void)" << endl; 
		cerr << "Error : an XMLException occurs during parsing " << fileName << endl;

		char* message = XMLString::transcode(toCatch.getMessage());
		cerr << "Error : " << message << endl;
		XMLString::release(&message);	

		exceptionOccurs = true;
	
	} catch(const DOMException & toCatch) {
	
		/////////////////////////////////////////////////////////////////////////////
		// Gestion des exceptions DOM.
		/////////////////////////////////////////////////////////////////////////////
		cerr << "Error : in gaussianBasisFunction::load(void)" << endl; 
		cerr << "Error : an DOMException occurs during parsing " << fileName << endl;
		
		char* message = XMLString::transcode(toCatch.msg);
		cerr << "Error : " << toCatch.getMessage() << endl;
		XMLString::release(&message); 	
	
		exceptionOccurs = true;
	} catch (...) {
		cerr << "Error : in gaussianBasisFunction::load(void)" << endl; 
		cerr << "Error : an exception occurs during parsing " << fileName << endl;
	 	exceptionOccurs = true;
	}
	
	//////////////////////////////////////////////////////////////////
	// Lorsque nous avons des erreurs on arrete le programe.
	//////////////////////////////////////////////////////////////////
	if(exceptionOccurs
	||
	errorReporter->getSawErrors())
	{
		exit(1);
	}
	
	////////////////////////////////////////////////////////////////////
	// On récupère ce qu'il faut pour le manippuler ensuite...
	////////////////////////////////////////////////////////////////////
	Document = Parser->getDocument();

	/////////////////////////////////////////////////////////////////////
	// On met le document root comme il faut.
	/////////////////////////////////////////////////////////////////////
	setRootElement(Document->getDocumentElement());
}

void xmlDocumentParser::setExternalNoNamespaceSchemaLocation(const string & schemaFileName) 
{
	//////////////////////////////////////////////////////////////////////////
	// Conversion de la chaine de charactère et on passe en argument.
	//////////////////////////////////////////////////////////////////////////
	XMLCh * xmlSchemaFileName = XMLString::transcode(schemaFileName.c_str());
	Parser->setExternalNoNamespaceSchemaLocation(xmlSchemaFileName);
	XMLString::release(&xmlSchemaFileName);

}

