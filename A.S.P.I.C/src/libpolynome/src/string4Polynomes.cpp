#include <dpoint.h>
#include <iostream>
#include <ipoint.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <string>
using namespace std;

////////////////////////////////////////////////////////////
// Méthode pour convertir un entier en chaine de charctère.
////////////////////////////////////////////////////////////
const char *  int2string(const int & i)
{
	static char charInt [64];
	sprintf(charInt,"%d",i);
	
	/////////////////////////////////////////////////
	// normalement on devrait catcher une exception.
	/////////////////////////////////////////////////

	return charInt;
}

////////////////////////////////////////////////////////////
// Méthode pour convertir un entier en chaine de charctère.
////////////////////////////////////////////////////////////
const char *  double2string(const double & value)
{
	static char charDouble [64];
	sprintf(charDouble,"%lf",value);
	
	/////////////////////////////////////////////////
	// normalement on devrait catcher une exception.
	/////////////////////////////////////////////////

	return charDouble;
}

////////////////////////////////////////////////////////////////
// Fonction pour érire un monome de la forme "(variable - centre)"
// - Lorsque le centre est nul on peut simplement écrire "varaible".
// - Lorsque le centre est non nul il faut prendre en compte le signe 
// poour éviter le (varialble - -1.000) par exemple.
////////////////////////////////////////////////////////////////
const string monome2string(const string & variable , const double & center)
{
	string base_string;
	
	if(center == 0) {
		base_string = variable;
		return base_string;
	}

	base_string = "(" + variable;
	
	if(center>0) {
		base_string += " - ";
	} else {
		base_string += " + ";
	}
	base_string += double2string(fabs(center));
	base_string += ")";

	return base_string;	
}

const string monome2string(const string & variable , const double & center ,const int & deg)
{
	if(deg == 0) {
		return "1";
	}
	
	string base_string = monome2string(variable,center);
	
	if(deg == 1) {
		return base_string;
	}
	
	base_string += "^";
	base_string += int2string(deg);
	
	return base_string;	
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////////////
const string monome3D2string(const string variables[3] , const ipoint<3> & degree , const dpoint<3> & center=dpoint<3>(0))
{
	if(degree == ipoint<3>(0)) {
		return "1";
	}

	string strMonome ="";

	for(int i=0 ; i < 3 ; i++) {
		if(degree[i] > 0) {
			strMonome += monome2string(variables[i],center[i],degree[i]);
		}
	}

	return strMonome;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
const string monome3D2string(const ipoint<3> & degree , const dpoint<3> & center=dpoint<3>(0))
{
	static string variables[3] = {"x","y","z"};
	return monome3D2string(variables,degree,center);
}

