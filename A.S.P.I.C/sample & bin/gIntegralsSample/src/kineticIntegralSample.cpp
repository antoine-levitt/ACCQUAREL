#include <monoElectronicIntegral.h>
#include <gaussianBasisFunction.h>
#include <xmlBasisFunctionDocumentParser.h>
#include <iostream>
#include "mapleDebugInterface.h"
using namespace std;




int main(int argc , char ** argv)
{
	double kinetic;
	
	gaussianBasisFunction phi_1 , phi_2;
	xmlBasisFunctionDocumentParser gbfParser;	
	monoElectronicIntegral kineticComputer; 
	
	////////////////////////////////////////////////////////////////////
	//  Validating command line arguments.
	////////////////////////////////////////////////////////////////////
	if(argc != 3) {
		cerr << "Error : in kinetic.exe" << endl;
		cerr << "Error : waiting for 2 arguments " << argc-1 << " where found." << endl;
		cerr << "Error : aborting process." << endl;
		return -1;
	}

	///////////////////////////////////////////////////////////////////
	//  Reading the first basis function.
	//////////////////////////////////////////////////////////////////	
	phi_1 = gbfParser.getBasisFunction4Document(argv[1]);	
	cout << phi_1 << endl;	
	
	///////////////////////////////////////////////////////////////////
	// Reading the second basis function.
	///////////////////////////////////////////////////////////////////
	phi_2 = gbfParser.getBasisFunction4Document(argv[2]);	
	cout << phi_2 << endl;	
	
	////////////////////////////////////////////////////////////////////
	// Calcul du recouvrement.
	////////////////////////////////////////////////////////////////////	
	kinetic =  kineticComputer.getKineticValue(phi_1,phi_2);
	cout << "Kinetic value is : " << setprecision(10) << kinetic << endl;	
	
	/////////////////////////////////////////////////////////////////////
	// Génération du fichier Mapple pour comparaison.
	//////////////////////////////////////////////////////////////////////
	mapleDebugInterface mKinetic;

	///////////////////////////////////////////////////////
	// 1 - génération du chemin vers le fichier.
	// et initialisation du flux.
	///////////////////////////////////////////////////////
	string filePath = getenv("ASPIC_ROOT");
	filePath += "/gaussian_integrals_sample/data/kinetic.txt";
	mKinetic.initOutFileStream(filePath);

	////////////////////////////////////////////////////////
	// 2 - On balance ici tout ce dont nous avons besoin.
	/////////////////////////////////////////////////////////
	string strBf1 , strBf2;

	mKinetic.getOutFileStream() << "> restart:" << endl;
	strBf1 = mKinetic.initMapleGaussianBasisFunction("1" , phi_1);
	strBf2 = mKinetic.initMapleGaussianBasisFunction("2" , phi_2);

	//////////////////////////////////////////////////////////
	// On fait le calcul.
	//////////////////////////////////////////////////////////
	mKinetic.getOutFileStream() << "> Kinetic := 0;" << endl;

	for(int i=1 ; i <= phi_1.getNbrOfContractions() ; i++) {
		for(int j=1 ; j <= phi_2.getNbrOfContractions() ; j++) {

		
			mKinetic.getOutFileStream() << "> Int(Int(Int(" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf1 << "[" << i << "](x,y,z) , x) *" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf2 << "[" << j << "](x,y,z) , x)" << endl;
			mKinetic.getOutFileStream() << "> ,x=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << "> ,y=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << ">, z=-infinity .. infinity);" << endl;

			mKinetic.getOutFileStream() << "> kinetic_temp := int(int(int(" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf1 << "[" << i << "](x,y,z) , x) *" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf2 << "[" << j << "](x,y,z) , x)" << endl;
			mKinetic.getOutFileStream() << "> ,x=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << "> ,y=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << ">, z=-infinity .. infinity);" << endl;

			mKinetic.getOutFileStream() << "> Kinetic := Kinetic + kinetic_temp;" << endl;

			mKinetic.getOutFileStream() << "> Int(Int(Int(" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf1 << "[" << i << "](x,y,z) , y) *" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf2 << "[" << j << "](x,y,z) , y)" << endl;
			mKinetic.getOutFileStream() << "> ,x=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << "> ,y=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << ">, z=-infinity .. infinity);" << endl;

			mKinetic.getOutFileStream() << "> kinetic_temp := int(int(int(" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf1 << "[" << i << "](x,y,z) , y) *" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf2 << "[" << j << "](x,y,z) , y)" << endl;
			mKinetic.getOutFileStream() << "> ,x=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << "> ,y=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << ">, z=-infinity .. infinity);" << endl;

			mKinetic.getOutFileStream() << "> Kinetic := Kinetic + kinetic_temp;" << endl;

			mKinetic.getOutFileStream() << "> Int(Int(Int(" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf1 << "[" << i << "](x,y,z) , z) *" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf2 << "[" << j << "](x,y,z) , z)" << endl;
			mKinetic.getOutFileStream() << "> ,x=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << "> ,y=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << ">, z=-infinity .. infinity);" << endl;

			mKinetic.getOutFileStream() << "> kinetic_temp := int(int(int(" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf1 << "[" << i << "](x,y,z) , z) *" << endl;
			mKinetic.getOutFileStream() << "> diff(" << strBf2 << "[" << j << "](x,y,z) , z)" << endl;
			mKinetic.getOutFileStream() << "> ,x=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << "> ,y=-infinity .. infinity)" << endl;
			mKinetic.getOutFileStream() << ">, z=-infinity .. infinity);" << endl;

			mKinetic.getOutFileStream() << "> Kinetic := Kinetic + kinetic_temp;" << endl;
		}
	}

	/////////////////////////////////////////////////////////////////////////
	// On affiche le résultat.
	/////////////////////////////////////////////////////////////////////////
	mKinetic.getOutFileStream() << "> evalf(Kinetic);" << endl;
	mKinetic.getOutFileStream() << "> Aspic_Kinetic := " << kinetic << ";" << endl;

	//////////////////////////////////////////////////////////
	// On arrete le flux.
	//////////////////////////////////////////////////////////
	mKinetic.terminateOutFileStream();

	return 0;
}

