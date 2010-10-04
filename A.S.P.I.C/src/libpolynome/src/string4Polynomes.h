#ifndef _STRING_4_POLYNOMES_
#define _STRING_4_POLYNOMES_

#include <ipoint.h>
#include <dpoint.h>

extern const char * int2string(const int & i);

extern const char * double2string(const double & i);

extern const string monome2string(const string & variable , const double & center);

extern const string monome2string(const string & variable , const double & center , const int & deg);

extern const string monome3D2string(const string variables[3] , const ipoint<3> & degree , const dpoint<3> & center=dpoint<3>(0));

extern const string monome3D2string(const ipoint<3> & degree , const dpoint<3> & center=dpoint<3>(0));

//extern bool testFilePath4Reading(const string & filPath);

#endif
