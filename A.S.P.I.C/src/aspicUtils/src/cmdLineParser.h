/* 
* The coniguration manager library of A.S.P.I.C. 
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
#ifndef _CMD_LINE_PARSER_
#define _CMD_LINE_PARSER_

#include <assert.h>
#include <map>
#include <string>
using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gestion des arguments ligne de commande.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
class cmdLineParser
{
	
private :

	
	string ProgramName;
	
	/**
	 * The map containing the arguments and options
	 * of the command line.
	 */
	map<string,string> Arguments;
	
protected:
	
	/**
	 * The copy method.
	 */
	void copy(const cmdLineParser & cmdLine);
	
public:
	
	/**
	 * The constructor.
	 */
	cmdLineParser(void);
	
	/**
	 * Constructor with cmd line arguments.
	 */
	cmdLineParser(int argc , char * argv[]);
	
	/**
	 * The copy constructor.
	 */
	cmdLineParser(const cmdLineParser & cmdLine);
	
	/**
		* The destructor.
	 */
	~cmdLineParser(void);

	/**
	 * An iterator to parse the map.
	 */
	map<string,string>::const_iterator begin(void) const;
	
	/**
		* The clear method.
	 */
	void clear(void);	
	
	/**
	 * Method to know if an object is empty.
	 */
	bool empty(void) const;
	
		/**
		* An iterator to parse the map.
		 */
	map<string,string>::const_iterator end(void) const;
	
	/**
	 * Method GET to find the option associeted with an argument.
	 */
	string getArgumentOption(const string & argName) const;
	
	/**
   * Method GET for all the arguments in the cmd line.
	 */
	const map<string,string> & getArguments(void) const;
	
	/**
	 * Method GET for the number of arguments.
	 */
	int getNbrOfArguments(void) const;
	
	/**
	 * Method GET for the name of the program.
	 */
	const string & getProgramName(void) const;
	
	/**
	 * Method to know if there is an argument with a name.
	 */
	bool hasArgument(const string & argName) const;
	
	/**
		* Method SET for a particular argument with option.
	 */
	map<string,string>::iterator setArgument(const string & argument , const string & option);
	
	/**
	 * Method that reads the command line.
	 */
	void setArguments(int argc , char * argv[]);

	/**
	 * Method SET for the arguments from a map.
	 */
	void setArguments(const map<string,string> & arguments);


	/**
	 * Method SET for the program name.
	 */
	void setProgramName(const string & programName);
};



#endif
