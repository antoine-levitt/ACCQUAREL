/*
 * Copyright 1999-2002,2004 The Apache Software Foundation.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// ---------------------------------------------------------------------------
//  Includes
// ---------------------------------------------------------------------------
#include <xercesc/sax/SAXParseException.hpp>
#include "xmlErrorReporter.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>
using namespace std;


void xmlErrorReporter::warning(const SAXParseException&)
{
    //
    // Ignore all warnings.
    //
}

void xmlErrorReporter::error(const SAXParseException& toCatch)
{
    fSawErrors = true;
		char * message = XMLString::transcode(toCatch.getSystemId()); 	
		cerr << "Error : in file \"" << message << "\"" << endl;
		XMLString::release(&message);

		cerr << "Error : line " << toCatch.getLineNumber() << " column " << toCatch.getColumnNumber() << endl;
    
		message = XMLString::transcode(toCatch.getMessage());
		cerr << "Error : " << message << endl;
		XMLString::release(&message);
}

void xmlErrorReporter::fatalError(const SAXParseException& toCatch)
{
    fSawErrors = true;
		char * message = XMLString::transcode(toCatch.getSystemId()); 	
		cerr << "Fatal Error : in file \"" << message << "\"" << endl;
		XMLString::release(&message);

		cerr << "Fatal Error : line " << toCatch.getLineNumber() << " column " << toCatch.getColumnNumber() << endl;
    
		message = XMLString::transcode(toCatch.getMessage());
		cerr << "Fatal Error : " << message << endl;
		XMLString::release(&message);
    
}

void xmlErrorReporter::resetErrors()
{
    fSawErrors = false;
}



