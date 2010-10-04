#include <aspicConfiguration.h>
#include <monoElectronicIntegral.h>
#include <gaussianBasisFunction.h>
#include <xmlBasisFunctionDocumentParser.h>
#include <iostream>
#include "mapleDebugInterface.h"
using namespace std;




int main(int argc , char ** argv)
{
	double overlap ;
	
	gaussianBasisFunction phi_1 , phi_2;
	xmlBasisFunctionDocumentParser gbfParser;	
	monoElectronicIntegral overlapComputer;
	
	////////////////////////////////////////////////////////////////////
	// Vérification des arguments passés en ligne de commande.
	////////////////////////////////////////////////////////////////////
	if(argc != 3) {
		cerr << "Error : in overlap.exe" << endl;
		cerr << "Error : waiting for 2 arguments " << argc-1 << " where found." << endl;
		cerr << "Error : aborting process." << endl;
		return -1;
	}

	///////////////////////////////////////////////////////////////////
	// Construction de la première fonction de base.
	//////////////////////////////////////////////////////////////////	
	phi_1 = gbfParser.getBasisFunction4Document(argv[1]);	
	cout << phi_1 << endl;	
	
	///////////////////////////////////////////////////////////////////
	// Construction de la seconde fonction de base.
	///////////////////////////////////////////////////////////////////
	phi_2 = gbfParser.getBasisFunction4Document(argv[2]);	
	cout << phi_2 << endl;	
	
	////////////////////////////////////////////////////////////////////
	// Calcul du recouvrement.
	////////////////////////////////////////////////////////////////////	
	overlap =  overlapComputer.getOverlapValue(phi_1,phi_2);
	cout << "La valeure du recouvrement est : " << setprecision(16) << overlap << endl;	
	//cout << "LA valeur de normalisation est de : " << 
	/////////////////////////////////////////////////////////////////////
	// Génération du fichier Mapple pour comparaison.
	//////////////////////////////////////////////////////////////////////
	mapleDebugInterface mOverlap;

	///////////////////////////////////////////////////////
	// 1 - génération du chemin vers le fichier.
	// et initialisation du flux.
	///////////////////////////////////////////////////////
	string filePath = aspicConfiguration::getAspicRoot();
	filePath += "/sample & bin/gIntegralsSample/data/overlap.txt";
	mOverlap.initOutFileStream(filePath);

	////////////////////////////////////////////////////////
	// 2 - On balance ici tout ce dont nous avons besoin.
	/////////////////////////////////////////////////////////
	string strBf1 , strBf2;

	mOverlap.getOutFileStream() << "> restart:" << endl;
	mOverlap.getOutFileStream() << "> Aspic_Overlap := " << setprecision(16) <<  overlap << ";" << endl;;
	strBf1 = mOverlap.initMapleGaussianBasisFunction("1" , phi_1);
	strBf2 = mOverlap.initMapleGaussianBasisFunction("2" , phi_2);

	//////////////////////////////////////////////////////////
	// On fait le calcul.
	//////////////////////////////////////////////////////////
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

	/////////////////////////////////////////////////////////////////////////
	// On affiche le résultat.
	/////////////////////////////////////////////////////////////////////////
	mOverlap.getOutFileStream() << "> evalf(Overlap);" << endl;
	mOverlap.getOutFileStream() << "> Aspic_Overlap;" << endl;
	//////////////////////////////////////////////////////////
	// On arrete le flux.
	//////////////////////////////////////////////////////////
	mOverlap.terminateOutFileStream();

	return 0;
}

