
#include <aspicConfiguration.h>
#include <gaussianBasisFunction.h>
#include <iostream>
#include <mapleDebugInterface.h>
#include <monoElectronicIntegral.h>
#include <xmlBasisFunctionDocumentParser.h>

using namespace std;


int main(int argc , char ** argv)
{
	double overlap ;
	gaussianBasisFunction phi_1 , phi_2;
	xmlBasisFunctionDocumentParser gbfParser;	
	monoElectronicIntegral overlapComputer;
	
	
	/////////////////////////////////////////////////////////////////////////////////////
	// this programs shall be working with 2 arguments : the path to XML files that contains
	// the basis functions. If the number of arguments is wrong, then the program is aborted
	// right here right now.
	/////////////////////////////////////////////////////////////////////////////////////
	if(argc != 3) {
		cerr << "Error : in overlap.exe" << endl;
		cerr << "Error : waiting for 2 arguments " << argc-1 << " where found." << endl;
		cerr << "Error : aborting process." << endl;
		return -1;
	}

	///////////////////////////////////////////////////////////////////
	// Reading the first basis function in the first argument.
	//////////////////////////////////////////////////////////////////	
	phi_1 = gbfParser.getBasisFunction4Document(argv[1]);	
	cout << phi_1 << endl;	
	
	///////////////////////////////////////////////////////////////////
	// Reading the second basis function in the second argument.
	///////////////////////////////////////////////////////////////////
	phi_2 = gbfParser.getBasisFunction4Document(argv[2]);	
	cout << phi_2 << endl;	
	
	////////////////////////////////////////////////////////////////////////////////////
	// Computation of the overlap value unsing the ASPIC monoElectronicIntegral class.
	////////////////////////////////////////////////////////////////////////////////////	
	overlap =  overlapComputer.getOverlapValue(phi_1,phi_2);
	cout << "La valeure du recouvrement est : " << setprecision(16) << overlap << endl;	

	
	/////////////////////////////////////////////////////////////////////
	// From now on the Maple file is written ...
	//////////////////////////////////////////////////////////////////////
	mapleDebugInterface mOverlap;

	
	string filePath = aspicConfiguration::getAspicRoot();
	filePath += "/sample & bin/gIntegralsSample/data/overlap.txt";
	mOverlap.initOutFileStream(filePath);

	string strBf1 , strBf2;

	mOverlap.getOutFileStream() << "> restart:" << endl;
	mOverlap.getOutFileStream() << "> Aspic_Overlap := " << setprecision(16) <<  overlap << ";" << endl;;
	strBf1 = mOverlap.initMapleGaussianBasisFunction("1" , phi_1);
	strBf2 = mOverlap.initMapleGaussianBasisFunction("2" , phi_2);

	mOverlap.getOutFileStream() << "> Overlap := 0;" << endl;

	for(int i=1 ; i <= phi_1.getNbrOfContractions() ; i++) {
		for(int j=1 ; j <= phi_2.getNbrOfContractions() ; j++) {

			mOverlap.getOutFileStream() << "> Int(Int(Int(" << endl;
			mOverlap.getOutFileStream() << "> " << strBf1 << "[" << i << "](x,y,z) *" << endl;
			mOverlap.getOutFileStream() << "> " << strBf2 << "[" << j << "](x,y,z)" << endl;
			mOverlap.getOutFileStream() << "> ,x=-infinity .. infinity)" << endl;
			mOverlap.getOutFileStream() << "> ,y=-infinity .. infinity)" << endl;
			mOverlap.getOutFileStream() << ">, z=-infinity .. infinity);" << endl;

			mOverlap.getOutFileStream() << "> overlap_temp := evalf[16](int(int(int(" << endl;
			mOverlap.getOutFileStream() << "> " << strBf1 << "[" << i << "](x,y,z) *" << endl;
			mOverlap.getOutFileStream() << "> " << strBf2 << "[" << j << "](x,y,z)" << endl;
			mOverlap.getOutFileStream() << "> ,x=-infinity .. infinity)" << endl;
			mOverlap.getOutFileStream() << "> ,y=-infinity .. infinity)" << endl;
			mOverlap.getOutFileStream() << ">, z=-infinity .. infinity));" << endl;

			mOverlap.getOutFileStream() << "> Overlap := Overlap + overlap_temp;" << endl;
		}
	}

	mOverlap.getOutFileStream() << "> evalf(Overlap);" << endl;
	mOverlap.getOutFileStream() << "> Aspic_Overlap;" << endl;
	mOverlap.terminateOutFileStream();

	// beast for the end : returning OK exit status.
	return 0;
}

