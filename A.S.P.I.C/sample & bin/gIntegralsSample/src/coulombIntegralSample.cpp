#include <aspicConfiguration.h>
#include <coulombIntegral.h>
#include <gaussianBasisFunction.h>
#include <hermiteProjector.h>
#include <hermiteIntegrationBuffer.h>
#include <iostream>
#include "mapleDebugInterface.h"
//#include <string_utils.h>
#include <xmlBasisFunctionDocumentParser.h>

using namespace std;


/**
 * On va tester la décomposition en polynomes de Hermite.
 */ 
void test_hermite_projector(const gaussianBasisFunction & phi_1 , const gaussianBasisFunction & phi_2)
{
	hermiteProjector hermite;

	// Ici on calcule la projection du produit des deux permiere contraction de nos fonction
	// de base sur une base de polyômes de Hermite adaptée.
	gaussianPolynome3D gp = phi_1.getGaussianPolynome3D(0) * phi_2.getGaussianPolynome3D(0);
	hermite.multiply(phi_1.getGaussianPolynome3D(0) , phi_2.getGaussianPolynome3D(0));

	// Création du fichier Maple pour aider au debbuguage.
	mapleDebugInterface mHermite;
	string filePath = aspicConfiguration::getAspicRoot();
	filePath += "/sample & bin/gIntegralsSample/data/hermiteProjector.txt";
	mHermite.initOutFileStream(filePath);
	mHermite.getOutFileStream() << "> restart;" << endl;

	// 1. on affiche le polynome produit basic.
	mHermite.initMaplePolynome("polynome_1" , phi_1.getGaussianPolynome3D(0));
	mHermite.initMaplePolynome("polynome_2" , phi_2.getGaussianPolynome3D(0));
	mHermite.initMaplePolynome("prodPolynome" , gp);

	mHermite.getOutFileStream() << " ###############################################" << endl;
	mHermite.getOutFileStream() << " # A cette ligne il faut trouver 0  !!!!!!!!!! #" << endl;
	mHermite.getOutFileStream() << " ###############################################" << endl;
	mHermite.getOutFileStream() << "> expand(prodPolynome - polynome_1 * polynome_2);" << endl;

	// 2. on affiche le polynome de hermite.
	
	ipoint<3> degree;
	dpoint<3> center = hermite.getCenter();
	double coefficient;

	mHermite.getOutFileStream() << "> hermitePoly_Exponent := " << hermite.getExponent() << ";" << endl;
	mHermite.getOutFileStream() << "> hermitePoly_Center_x := " << hermite.getCenter()[0] << ";" << endl;
	mHermite.getOutFileStream() << "> hermitePoly_Center_y := " << hermite.getCenter()[1] << ";" << endl;
	mHermite.getOutFileStream() << "> hermitePoly_Center_z := " << hermite.getCenter()[2] << ";" << endl;
	
	mHermite.getOutFileStream() << "> hermitePoly := " << endl; 
	for(degree = hermite.getHermiteFirst() ; degree != hermite.getHermiteEnd() ; degree = hermite.getHermiteNext(degree)) {
		coefficient = hermite.getHermiteCoeffcient(degree);
		mHermite.getOutFileStream() << "> " << (coefficient > 0 ? "+" : "") << coefficient << " *";
		mHermite.getOutFileStream() << "sqrt(hermitePoly_Exponent)^" << degree.norme_1() << " *" ; 
		mHermite.getOutFileStream() << "HermiteH(" << degree[0] <<" , sqrt(hermitePoly_Exponent)*(x-hermitePoly_Center_x)) *" ; 
		mHermite.getOutFileStream() << "HermiteH(" << degree[1] <<" , sqrt(hermitePoly_Exponent)*(y-hermitePoly_Center_y)) *" ; 
		mHermite.getOutFileStream() << "HermiteH(" << degree[2] <<" , sqrt(hermitePoly_Exponent)*(z-hermitePoly_Center_z)) " << endl ; 
	}
	mHermite.getOutFileStream() << ">;" << endl; 

	mHermite.terminateOutFileStream();
}


/**
 * Cette fonction permet de vérifier que la classe 
 * hermiteBasisProjector calcule bien les polynomes
 * que nous cherchons.
 */
void testHermiteBasisPolynome(void)
{
	hermiteBasisPolynome hermiteBasisPoly;
	mapleDebugInterface mHermiteBasisPoly;
	int degree , degree_max = 40;

	hermiteBasisPoly.setCenter(3);
	hermiteBasisPoly.setExponent(1/3.);
	hermiteBasisPoly.setDegree(0);
	
	string filePath = getenv("ASPIC_ROOT");
	filePath += "/gaussian_integrals_sample/data/hermiteBasisPolynome.txt";
	mHermiteBasisPoly.initOutFileStream(filePath);
	mHermiteBasisPoly.getOutFileStream() << "> restart;" << endl;
	
	mHermiteBasisPoly.getOutFileStream() << "> R := " << hermiteBasisPoly.getCenter() << ";" << endl;
	mHermiteBasisPoly.getOutFileStream() << "> v := sqrt(" << hermiteBasisPoly.getExponent() << ");" << endl;

	for(degree = 0 ; degree < (degree_max+1) ; degree++ ) {
		string polyBaseName = "Hermite";
		polyBaseName += int2string(degree);
		mHermiteBasisPoly.initMaplePolynome(polyBaseName , hermiteBasisPoly.getPolynome());
		mHermiteBasisPoly.getOutFileStream() << "> v^" << degree << " * HermiteH(" << degree <<" , v*R*x):" << endl;
		mHermiteBasisPoly.getOutFileStream() << "> expand(% - " << polyBaseName<< ");" << endl;
	
		if(degree < degree_max ) {
			hermiteBasisPoly.addDegree();
		}
	}

	mHermiteBasisPoly.terminateOutFileStream();
}


/**
 * Test du buffer d'intégration.
 */
void testHermiteIntegrationBuffer(const double & v2 , const dpoint<3> R , const ipoint<3> & degreeMax)
{
	hermiteIntegrationBuffer buffer;
	mapleDebugInterface mBuffer;

	// On construit l'interface.
	string filePath = getenv("ASPIC_ROOT");
	filePath += "/gaussian_integrals_sample/data/hermiteIntegrationBuffer.txt";
	mBuffer.initOutFileStream(filePath);
	mBuffer.getOutFileStream() << "> restart;" << endl;

	// On construit les projection 
	buffer.setCenter(R);
	buffer.setExponent(v2);
	buffer.setIntegrationBufferSize(degreeMax);
	
	// Initilisation des variables Maple.
	mBuffer.getOutFileStream() <<"> v := sqrt(" << v2 << ");" << endl;
	mBuffer.getOutFileStream() <<"> Rx := " << R[0] << ";" << endl;
	mBuffer.getOutFileStream() <<"> Ry := " << R[1] << ";" << endl;
	mBuffer.getOutFileStream() <<"> Rz := " << R[2] << ";" << endl;


	for( int rx = 0 ; rx < buffer.getIntegrationBufferSize(0) ; rx++) {
		for( int ry = 0 ; ry < buffer.getIntegrationBufferSize(1) ; ry++) {
			for( int rz = 0 ; rz < buffer.getIntegrationBufferSize(2) ; rz++) {
				
				mBuffer.getOutFileStream() << "> Int (" << endl;
				mBuffer.getOutFileStream() << "> (v*u)^(" << rx +ry +rz << ") * " << endl;
				mBuffer.getOutFileStream() << "> HermiteH(" << rx << " , v*Rx*u) * " << endl;
				mBuffer.getOutFileStream() << "> HermiteH(" << ry << " , v*Ry*u) * " << endl;
				mBuffer.getOutFileStream() << "> HermiteH(" << rz << " , v*Rz*u) * " << endl;
				mBuffer.getOutFileStream() << "> exp( - v^2 * (Rx^2 + Ry^2 + Rz^2) * u^2)" << endl;
				mBuffer.getOutFileStream() << "> , u=0..1);" <<endl;
				
				mBuffer.getOutFileStream() << "> maple_value := evalf(int (" << endl;
				mBuffer.getOutFileStream() << "> (v*u)^(" << rx +ry +rz << ") * " << endl;
				mBuffer.getOutFileStream() << "> HermiteH(" << rx << " , v*Rx*u) * " << endl;
				mBuffer.getOutFileStream() << "> HermiteH(" << ry << " , v*Ry*u) * " << endl;
				mBuffer.getOutFileStream() << "> HermiteH(" << rz << " , v*Rz*u) * " << endl;
				mBuffer.getOutFileStream() << "> exp( - v^2 * (Rx^2 + Ry^2 + Rz^2) * u^2)" << endl;
				mBuffer.getOutFileStream() << "> , u=0..1));" <<endl;
				
				mBuffer.getOutFileStream() << "> aspic_value := " << setprecision(15) << buffer.getIntegrationValue(rx,ry,rz) << ";" << endl;
				mBuffer.getOutFileStream() << "> error_value := maple_value - aspic_value;" << endl;

				
			}
		}
	}
}


void testCoulomb(const gaussianBasisFunction & phi_1 , const gaussianBasisFunction & phi_2 , const gaussianBasisFunction & phi_3 , const gaussianBasisFunction & phi_4)
{
	coulombIntegral qe;
	double aspic_coulomb = qe.getCoulombValue(phi_1,phi_2,phi_3,phi_4);


	///////////// ici on écrit un programe mapple pour faire le calcul ///////////////////////////
	//hermiteIntegrationBuffer buffer;
	hermiteProjector hermite_12 , hermite_34;
	mapleDebugInterface mBuffer;

	// On construit l'interface.
	string filePath = aspicConfiguration::getAspicRoot();
	filePath += "/sample & bin/gIntegralsSample/data/coulomb.txt";
	mBuffer.initOutFileStream(filePath);
	mBuffer.getOutFileStream() << "> restart;" << endl;
	mBuffer.getOutFileStream() << "> maple_coulomb := 0;" << endl;
	mBuffer.getOutFileStream() << "> aspic_coulomb := " << setprecision(16) << aspic_coulomb << ";" << endl;
	

	int i , j , k , l , Ka , Kb , Kc, Kd;
	gaussianPolynome3D gp_a, gp_b , gp_c , gp_d;

	Ka = phi_1.getNbrOfContractions();
	Kb = phi_2.getNbrOfContractions();
	Kc = phi_3.getNbrOfContractions();
	Kd = phi_4.getNbrOfContractions();


	// On boucle sur le nombre de contractions de la première
	// fonction de base.
	for(i=0 ; i < Ka ; i++) {
		for(j=0 ; j < Kb ; j++) {
			for(k=0 ; k < Kc ; k++) {				
				for(l=0 ; l < Kd ; l++) {

					// On construit les projection 
					hermite_12.multiply(phi_1.getGaussianPolynome3D(i) , phi_2.getGaussianPolynome3D(j));
					hermite_34.multiply(phi_3.getGaussianPolynome3D(k) , phi_4.getGaussianPolynome3D(l));
	
					mBuffer.getOutFileStream() <<"> Cab := "  <<  hermite_12.getCoefficient()  << ";" << endl;
					mBuffer.getOutFileStream() <<"> alpha := "  <<  hermite_12.getExponent()  << ";" << endl;
					mBuffer.getOutFileStream() <<"> Cabx := " <<  hermite_12.getCenter()[0] << ";" << endl;
					mBuffer.getOutFileStream() <<"> Caby := " <<  hermite_12.getCenter()[1] << ";" << endl;
					mBuffer.getOutFileStream() <<"> Cabz := " <<  hermite_12.getCenter()[2] << ";" << endl;
	
					mBuffer.getOutFileStream() <<"> Ccd	 := "  <<  hermite_34.getCoefficient()  << ";" << endl;
					mBuffer.getOutFileStream() <<"> beta := "  <<  hermite_34.getExponent()  << ";" << endl;
					mBuffer.getOutFileStream() <<"> Ccdx := " <<  hermite_34.getCenter()[0] << ";" << endl;
					mBuffer.getOutFileStream() <<"> Ccdy := " <<  hermite_34.getCenter()[1] << ";" << endl;
					mBuffer.getOutFileStream() <<"> Ccdz := " <<  hermite_34.getCenter()[2] << ";" << endl;

					// Initilisation des variables Maple.
					mBuffer.getOutFileStream() <<"> C := 2*Pi^(5/2) / (alpha * beta * sqrt(alpha+beta)) * Cab * Ccd;" << endl;
					mBuffer.getOutFileStream() <<"> v := sqrt(alpha * beta / (alpha+beta));" << endl;
					mBuffer.getOutFileStream() <<"> Rx := -Cabx + Ccdx;" << endl;
					mBuffer.getOutFileStream() <<"> Ry := -Caby + Ccdy;" << endl;
					mBuffer.getOutFileStream() <<"> Rz := -Cabz + Ccdz;" << endl;
					mBuffer.getOutFileStream() <<"> T := v^2 * (Rx^2 + Ry^2 +Rz^2);" << endl;

					
					mBuffer.getOutFileStream() << "> maple_partial_coulomb := 0;" << endl;
					for(ipoint<3> i = hermite_12.getHermiteFirst() ; i != hermite_12.getHermiteEnd() ; i = hermite_12.getHermiteNext(i)) {
						for(ipoint<3> j = hermite_34.getHermiteFirst() ; j != hermite_34.getHermiteEnd() ; j = hermite_34.getHermiteNext(j)) {
							mBuffer.getOutFileStream() << "> maple_partial_coulomb := maple_partial_coulomb +" << endl;
							mBuffer.getOutFileStream() <<"> (" << setprecision(16) <<  hermite_12.getHermiteCoeffcient(i) * hermite_34.getHermiteCoeffcient(j) << ")" << endl;
							mBuffer.getOutFileStream() <<"> *(-1)^(" << j.norme_1() << ")* Int(" << endl;
							mBuffer.getOutFileStream() <<"> (v*u)^(" << i[0] << "+" << j[0] << ")*HermiteH(" << j[0] << "+" << i[0] << ",Rx * v* u)*" << endl;
							mBuffer.getOutFileStream() <<"> (v*u)^(" << i[1] << "+" << j[1] << ")*HermiteH(" << j[1] << "+" << i[1] << ",Ry * v* u)*" << endl;
							mBuffer.getOutFileStream() <<"> (v*u)^(" << i[2] << "+" << j[2] << ")*HermiteH(" << j[2] << "+" << i[2] << ",Rz * v* u)*" << endl;
							mBuffer.getOutFileStream() <<"> exp(-T*u^2),u=0..1)" << endl;
							mBuffer.getOutFileStream() << ">;" << endl;

						}
					}
					mBuffer.getOutFileStream() << "> maple_coulomb := maple_coulomb + C*maple_partial_coulomb;" << endl;
				}
			}
		}
	}
	
	mBuffer.getOutFileStream() << "> evalf[16](maple_coulomb);" << endl;
	mBuffer.getOutFileStream() << "> evalf[16](maple_coulomb - aspic_coulomb);" << endl;

}

int main(int argc , char ** argv)
{
	gaussianBasisFunction phi_1 , phi_2 , phi_3 , phi_4;
	xmlBasisFunctionDocumentParser gbfParser;	

	////////////////////////////////////////////////////////////////////
	// Vérification des arguments passés en ligne de commande.
	////////////////////////////////////////////////////////////////////
	if(argc != 5) {
		cerr << "Error : in overlap.exe" << endl;
		cerr << "Error : waiting for 4 arguments " << argc-1 << " where found." << endl;
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
	
	///////////////////////////////////////////////////////////////////
	// Construction de la seconde fonction de base.
	///////////////////////////////////////////////////////////////////
	phi_3 = gbfParser.getBasisFunction4Document(argv[3]);	
	cout << phi_3 << endl;	
	
	///////////////////////////////////////////////////////////////////
	// Construction de la seconde fonction de base.
	///////////////////////////////////////////////////////////////////
	phi_4 = gbfParser.getBasisFunction4Document(argv[4]);	
	cout << phi_4 << endl;	
		
	testCoulomb(phi_1,phi_2,phi_3,phi_4);
	//test_hermite_projector(phi_1,phi_2);
	//test_hermite_projector(phi_3,phi_4);

	//testHermiteIntegrationBuffer(2,dpoint<3>(0),ipoint<3>(10));
	return 0;
}

