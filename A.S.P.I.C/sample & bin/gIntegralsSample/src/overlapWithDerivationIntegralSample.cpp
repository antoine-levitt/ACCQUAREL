#include <cmdLineParser.h>
#include <gaussianBasisFunction.h>
#include <ipoint.h>
#include <iostream>
#include "mapleDebugInterface.h"
#include <monoElectronicIntegral.h>
#include <xmlBasisFunctionDocumentParser.h>
#include <xmlIPoint3DocumentParser.h>
using namespace std;



//////////////////////////////////////////////////////////////////////////////////////////////////
// Function that prints the options available for this program.
//
//////////////////////////////////////////////////////////////////////////////////////////////////
void printUsage(ostream & out) 
{
	out << endl;
	out << "Usage for overlapWithDeraivationIntegralDebugSample.exe" << endl;
	out << endl;
	out << "--phi1 <basisFunction.xml> : the first basis function" << endl;
	out << "--phi2 <basisFunction.xml> : the second basis function" << endl;
	out << "--deg1 <degree.xml>        : the derivation degree for the first basis function" << endl;
	out << "--deg2 <degree.xml>        : the derivation degree for the second basis function" << endl;
	out << "--output <filename.txt>    : path to the maple file that shall be written." << endl;
	out << "--help                     : print this message." << endl;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// Function that really computes the overlap with derivatives and generates the 
// maple text file.
//////////////////////////////////////////////////////////////////////////////////////////////////
void mapleDebugOverlap(const gaussianBasisFunction & phi1 ,
											 const ipoint<3> & deg1 ,
											 const gaussianBasisFunction & phi2 ,
											 const ipoint<3> & deg2,
											 const string & mapleFileName)
{
	////////////////////////////////////////////////////////////////////
	// Computation of the overlap value with derivation using the ASPIC
	// classes.
	////////////////////////////////////////////////////////////////////
	
	monoElectronicIntegral overlap;
	double overlapValue = overlap.getOverlapValue(phi1,deg1,phi2,deg2);
	
	//////////////////////////////////////////////////////////////////////
	// Generating maple text file for debugging.
	//////////////////////////////////////////////////////////////////////
	mapleDebugInterface mapleDebug;
	
	mapleDebug.initOutFileStream(mapleFileName);
	mapleDebug.getOutFileStream() << "> AspicOverlap := " << overlapValue << ";" << endl;
	
	string strBf1 = mapleDebug.initMapleGaussianBasisFunction("Phi1",phi1);
	string strBf2 = mapleDebug.initMapleGaussianBasisFunction("Phi2",phi2);
	
	mapleDebug.getOutFileStream() << "> Derivate1_x := " << deg1[0] << ";" << endl;
	mapleDebug.getOutFileStream() << "> Derivate1_y := " << deg1[1] << ";" << endl;
	mapleDebug.getOutFileStream() << "> Derivate1_z := " << deg1[2] << ";" << endl;
	mapleDebug.getOutFileStream() << "> Derivate2_x := " << deg2[0] << ";" << endl;
	mapleDebug.getOutFileStream() << "> Derivate2_y := " << deg2[1] << ";" << endl;
	mapleDebug.getOutFileStream() << "> Derivate2_z := " << deg2[2] << ";" << endl;

		mapleDebug.getOutFileStream() << "> Overlap := 0;" << endl;
		
		for(int i=1 ; i <= phi1.getNbrOfContractions() ; i++) {
			for(int j=1 ; j <= phi2.getNbrOfContractions() ; j++) {
				
				mapleDebug.getOutFileStream() << "> Int(Int(Int(" << endl;
				mapleDebug.getOutFileStream() << "> Diff(" << endl;
				mapleDebug.getOutFileStream() << "> " << strBf1 << "[" << i << "](x,y,z) , " << endl;
				mapleDebug.getOutFileStream() << "> x$Derivate1_x , y$Derivate1_y , z$Derivate1_z) *" << endl;
				mapleDebug.getOutFileStream() << "> Diff(" << endl;
				mapleDebug.getOutFileStream() << "> " << strBf2 << "[" << j << "](x,y,z) , " << endl;
				mapleDebug.getOutFileStream() << "> x$Derivate2_x , y$Derivate2_y , z$Derivate2_z) *" << endl;
				mapleDebug.getOutFileStream() << "> ,x=-infinity .. infinity)" << endl;
				mapleDebug.getOutFileStream() << "> ,y=-infinity .. infinity)" << endl;
				mapleDebug.getOutFileStream() << ">, z=-infinity .. infinity);" << endl;
				
				mapleDebug.getOutFileStream() << "> overlap_temp := evalf[16](int(int(int(" << endl;
				mapleDebug.getOutFileStream() << "> diff(" << endl;
				mapleDebug.getOutFileStream() << "> " << strBf1 << "[" << i << "](x,y,z) , " << endl;
				mapleDebug.getOutFileStream() << "> x$Derivate1_x , y$Derivate1_y , z$Derivate1_z) *" << endl;
				mapleDebug.getOutFileStream() << "> diff(" << endl;
				mapleDebug.getOutFileStream() << "> " << strBf2 << "[" << j << "](x,y,z) , " << endl;
				mapleDebug.getOutFileStream() << "> x$Derivate2_x , y$Derivate2_y , z$Derivate2_z) *" << endl;
				mapleDebug.getOutFileStream() << "> ,x=-infinity .. infinity)" << endl;
				mapleDebug.getOutFileStream() << "> ,y=-infinity .. infinity)" << endl;
				mapleDebug.getOutFileStream() << ">, z=-infinity .. infinity);" << endl;
				
				mapleDebug.getOutFileStream() << "> Overlap := Overlap + overlap_temp;" << endl;
			}
		}
		

}



int main(int argc , char ** argv)
{
	
	cmdLineParser cmdLine(argc,argv);
	gaussianBasisFunction phi1 , phi2;
	ipoint<3> deg1 , deg2;
	xmlBasisFunctionDocumentParser basisParser;
	xmlIPoint3DocumentParser degParser;
	string fileName;
	
	////////////////////////////////////////////////////////////////////////////////////////////
	// Help ?
	////////////////////////////////////////////////////////////////////////////////////////////
	if(cmdLine.hasArgument("help")) {
		printUsage(cout);
		return 0;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////
	// The  first gaussian basis function :
	////////////////////////////////////////////////////////////////////////////////////////////
	if(cmdLine.hasArgument("phi1") && !(fileName = cmdLine.getArgumentOption("phi1")).empty()) {
		phi1 = basisParser.getBasisFunction4Document(fileName);
		basisParser.close();
	} else {
		cerr << "---------------------------------------------------------------------------" << endl;
		cerr << "Error : in overlapWithDerivativeSample.exe" << endl;
		cerr << "Error : argument \"phi1\" not found." << endl;
		cerr << "Error : Option for the overlapWithDerivativeSample.exe are displayed below." << endl;
		cerr << "---------------------------------------------------------------------------" << endl;
		printUsage(cerr);
		cerr << "---------------------------------------------------------------------------" << endl;
		exit(1);
	}
	
	clog << "---------------------------------------------------------------------------" << endl;
	clog << " The first basis function " << endl;
	clog << "---------------------------------------------------------------------------" << endl;
	clog << phi1 << endl;
	clog << "---------------------------------------------------------------------------" << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////
	// The  second gaussian basis function :
	////////////////////////////////////////////////////////////////////////////////////////////
	if(cmdLine.hasArgument("phi2") && !(fileName = cmdLine.getArgumentOption("phi2")).empty()) {
		phi2 = basisParser.getBasisFunction4Document(fileName);
		basisParser.close();
	} else {
		cerr << "---------------------------------------------------------------------------" << endl;
		cerr << "Error : in overlapWithDerivativeSample.exe" << endl;
		cerr << "Error : argument \"phi2\" not found." << endl;
		cerr << "Error : Option for the overlapWithDerivativeSample.exe are displayed below." << endl;
		cerr << "---------------------------------------------------------------------------" << endl;
		printUsage(cerr);
		cerr << "---------------------------------------------------------------------------" << endl;
		exit(1);
	}
		
	clog << "---------------------------------------------------------------------------" << endl;
	clog << " The second basis function " << endl;
	clog << "---------------------------------------------------------------------------" << endl;
	clog << phi2 << endl;
	clog << "---------------------------------------------------------------------------" << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////
	// The first derivation degree :
	////////////////////////////////////////////////////////////////////////////////////////////
	if(cmdLine.hasArgument("deg1") && !(fileName = cmdLine.getArgumentOption("deg1")).empty()) {
		deg1 = degParser.getIPoint34Document(fileName);
		degParser.close();
	} else {
		cerr << "---------------------------------------------------------------------------" << endl;
		cerr << "Error : in overlapWithDerivativeSample.exe" << endl;
		cerr << "Error : argument \"deg1\" not found." << endl;
		cerr << "Error : Option for the overlapWithDerivativeSample.exe are displayed below." << endl;
		cerr << "---------------------------------------------------------------------------" << endl;
		printUsage(cerr);
		cerr << "---------------------------------------------------------------------------" << endl;
		exit(1);
	}
	
	clog << "---------------------------------------------------------------------------" << endl;
	clog << " The first derivation degree " << endl;
	clog << "---------------------------------------------------------------------------" << endl;
	clog << deg1 << endl;
	clog << "---------------------------------------------------------------------------" << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////
	// The first derivation degree :
	////////////////////////////////////////////////////////////////////////////////////////////
	if(cmdLine.hasArgument("deg2") && !(fileName = cmdLine.getArgumentOption("deg2")).empty()) {
		deg2 = degParser.getIPoint34Document(fileName);
		degParser.close();
	} else {
		cerr << "---------------------------------------------------------------------------" << endl;
		cerr << "Error : in overlapWithDerivativeSample.exe" << endl;
		cerr << "Error : argument \"deg2\" not found." << endl;
		cerr << "Error : Option for the overlapWithDerivativeSample.exe are displayed below." << endl;
		cerr << "---------------------------------------------------------------------------" << endl;
		printUsage(cerr);
		cerr << "---------------------------------------------------------------------------" << endl;
		exit(1);
	}
	
	
	////////////////////////////////////////////////////////////////////////////////////////////
	// The first derivation degree :
	////////////////////////////////////////////////////////////////////////////////////////////
	if(cmdLine.hasArgument("output") && !(fileName = cmdLine.getArgumentOption("output")).empty()) {
		;
	} else {
		cerr << "---------------------------------------------------------------------------" << endl;
		cerr << "Error : in overlapWithDerivativeSample.exe" << endl;
		cerr << "Error : argument \"output\" not found." << endl;
		cerr << "Error : Option for the overlapWithDerivativeSample.exe are displayed below." << endl;
		cerr << "---------------------------------------------------------------------------" << endl;
		printUsage(cerr);
		cerr << "---------------------------------------------------------------------------" << endl;
		exit(1); 
	}
	
	
	mapleDebugOverlap(phi1 , deg1 , phi2 , deg2, fileName);
	
	return 0;
}

