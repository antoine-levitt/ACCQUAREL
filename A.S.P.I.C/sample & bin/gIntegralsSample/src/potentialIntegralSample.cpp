#include <aspicConfiguration.h>
#include <gaussianBasisFunction.h>
#include <hermiteProjector.h>
#include <hermiteIntegrationBuffer.h>
#include <iostream>
#include "mapleDebugInterface.h"
#include <potentialIntegral.h>
//#include <string_utils.h>
#include <xmlDPoint3DocumentParser.h>
#include <xmlBasisFunctionDocumentParser.h>

using namespace std;

void testPotential(const gaussianBasisFunction & phi_1 , const gaussianBasisFunction & phi_2 , const dpoint<3> & center)
{

	potentialIntegral potential;
	double aspic_potential = potential.getPotentialValue(phi_1,phi_2,center,atom::ATOMIC_UNIT);

	//hermiteIntegrationBuffer buffer;
	hermiteProjector hermite_12;

	mapleDebugInterface mBuffer;

	// On construit l'interface.
	string filePath = aspicConfiguration::getAspicRoot();
	filePath += "/gaussian_integrals_sample/data/potential.txt";
	mBuffer.initOutFileStream(filePath);
	mBuffer.getOutFileStream() << "> restart;" << endl;
	mBuffer.getOutFileStream() << "> maple_potential := 0;" << endl;
	mBuffer.getOutFileStream() << "> aspic_potential := " << setprecision(16) << aspic_potential << ";" << endl;
	

	int i , j , Ka , Kb ;
	gaussianPolynome3D gp_a, gp_b;

	Ka = phi_1.getNbrOfContractions();
	Kb = phi_2.getNbrOfContractions();


	mBuffer.getOutFileStream() <<"> Center_x := " <<  center[0] << ";" << endl;
	mBuffer.getOutFileStream() <<"> Center_y := " <<  center[1] << ";" << endl;
	mBuffer.getOutFileStream() <<"> Center_z := " <<  center[2] << ";" << endl;

	// On boucle sur le nombre de contractions de la première
	// fonction de base.
	for(i=0 ; i < Ka ; i++) {
		for(j=0 ; j < Kb ; j++) {

					// On construit les projection 
					hermite_12.multiply(phi_1.getGaussianPolynome3D(i) , phi_2.getGaussianPolynome3D(j));
	
					mBuffer.getOutFileStream() <<"> Cab := "  <<  hermite_12.getCoefficient()  << ";" << endl;
					mBuffer.getOutFileStream() <<"> alpha := "  <<  hermite_12.getExponent()  << ";" << endl;
					mBuffer.getOutFileStream() <<"> Cabx := " <<  hermite_12.getCenter()[0] << ";" << endl;
					mBuffer.getOutFileStream() <<"> Caby := " <<  hermite_12.getCenter()[1] << ";" << endl;
					mBuffer.getOutFileStream() <<"> Cabz := " <<  hermite_12.getCenter()[2] << ";" << endl;
	

					// Initilisation des variables Maple.
					mBuffer.getOutFileStream() <<"> C := 2*Pi /alpha * Cab;" << endl;
					mBuffer.getOutFileStream() <<"> v := sqrt(alpha);" << endl;
					mBuffer.getOutFileStream() <<"> Rx := -Cabx + Center_x;" << endl;
					mBuffer.getOutFileStream() <<"> Ry := -Caby + Center_y;" << endl;
					mBuffer.getOutFileStream() <<"> Rz := -Cabz + Center_z;" << endl;
					mBuffer.getOutFileStream() <<"> T := v^2 * (Rx^2 + Ry^2 +Rz^2);" << endl;

					
					mBuffer.getOutFileStream() << "> maple_partial_potential := 0;" << endl;
					for(ipoint<3> i = hermite_12.getHermiteFirst() ; i != hermite_12.getHermiteEnd() ; i = hermite_12.getHermiteNext(i)) {
							mBuffer.getOutFileStream() << "> maple_partial_potential := maple_partial_potential +" << endl;
							mBuffer.getOutFileStream() <<"> (" << setprecision(16) <<  hermite_12.getHermiteCoeffcient(i) << ")" << endl;
							mBuffer.getOutFileStream() <<"> * Int(" << endl;
							mBuffer.getOutFileStream() <<"> (v*u)^(" << i[0] << ")*HermiteH(" << i[0] << ",Rx * v* u)*" << endl;
							mBuffer.getOutFileStream() <<"> (v*u)^(" << i[1] << ")*HermiteH(" << i[1] << ",Ry * v* u)*" << endl;
							mBuffer.getOutFileStream() <<"> (v*u)^(" << i[2] << ")*HermiteH(" << i[2] << ",Rz * v* u)*" << endl;
							mBuffer.getOutFileStream() <<"> exp(-T*u^2),u=0..1)" << endl;
							mBuffer.getOutFileStream() << ">;" << endl;
					}

					mBuffer.getOutFileStream() << "> maple_potential := maple_potential + C*maple_partial_potential;" << endl;
				}
			}
	
	mBuffer.getOutFileStream() << "> evalf[16](maple_potential);" << endl;
	mBuffer.getOutFileStream() << "> evalf[16](maple_potential - aspic_potential);" << endl;

}

int main(int argc , char ** argv)
{
	gaussianBasisFunction phi_1 , phi_2;
	xmlBasisFunctionDocumentParser gbfParser;	
	xmlDPoint3DocumentParser dPoint3Parser;
	dpoint<3> center;

	////////////////////////////////////////////////////////////////////
	// Validating the command line
	////////////////////////////////////////////////////////////////////
	if(argc != 4) {
		cerr << "Error : in overlap.exe" << endl;
		cerr << "Error : waiting for 3 arguments " << argc-1 << " where found." << endl;
		cerr << "Error : aborting process." << endl;
		return -1;
	}

	///////////////////////////////////////////////////////////////////
	// The first basis function.
	//////////////////////////////////////////////////////////////////	
	phi_1 = gbfParser.getBasisFunction4Document(argv[1]);	
	cout << phi_1 << endl;	
	
	///////////////////////////////////////////////////////////////////
	// The second basis function.
	///////////////////////////////////////////////////////////////////
	phi_2 = gbfParser.getBasisFunction4Document(argv[2]);	
	cout << phi_2 << endl;	
	
	center = dPoint3Parser.getDPoint34Document(argv[3]);

	testPotential(phi_1,phi_2,center);
	return 0;
}

