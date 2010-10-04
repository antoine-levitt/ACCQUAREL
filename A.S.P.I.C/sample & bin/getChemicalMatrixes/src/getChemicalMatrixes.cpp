#include <monoElectronicIntegral.h>
#include "cmdLineParser.h"
#include <coulombIntegral.h>
#include <fstream>
#include <gaussianBasisFunction.h>
#include <moleculeMapper.h>
#include <xmlMoleculeDocumentParser.h>
#include <potentialIntegral.h>
//#include <string_utils.h>
#include <time.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function : humanPrintOverlap.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void humanPrintOverlap( const molecule & mol , ostream & outFileStream , int width = 15 , int precision = 5) 
{
	int nbrOfBasisFunctions , i, j;
	moleculeMapper map;
	monoElectronicIntegral overlap;
	gaussianBasisFunction phi_a , phi_b;
	double value;
	time_t start,end;

	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();

	time(&start);
	outFileStream << "---------------------------------------------------------------" << endl; 	
	outFileStream << "-                      Begin Overlap                          -" << endl;
	outFileStream << "---------------------------------------------------------------" << endl; 	
	for(i=0 ; i < nbrOfBasisFunctions ; i++) {		
		phi_a = map.getBasisFunction4Global(i);
		
		for(j=0 ; j < nbrOfBasisFunctions ; j++) {
			
			if(j<i) {
				outFileStream << setw(width) << setprecision(precision) <<"x";
				
			} else {
				phi_b = map.getBasisFunction4Global(j);

				value = overlap.getOverlapValue(phi_a,phi_b);
				outFileStream << setw(width) << setprecision(precision) << value;
			}	
		}
		outFileStream << endl;
	}
	outFileStream << "---------------------------------------------------------------" << endl; 	
	outFileStream << "-                        End Overlap                          -" << endl;
	outFileStream << "---------------------------------------------------------------" << endl; 	
	time(&end);
	outFileStream << "- Computation time : " << difftime(end,start) << endl;
	outFileStream << "---------------------------------------------------------------" << endl << endl << endl; 	

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function : humanPrintHamiltonian.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void humanPrintKinetic( const molecule & mol , ostream & outFileStream , int width = 15 , int precision = 5) 
{

	int nbrOfBasisFunctions , i, j;
	moleculeMapper map;
	monoElectronicIntegral overlap;
	gaussianBasisFunction phi_a , phi_b;
	double value;
	time_t start,end;

	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();

	time(&start);
	outFileStream << "---------------------------------------------------------------" << endl; 	
	outFileStream << "-                      Begin Kinetic                          -" << endl;
	outFileStream << "---------------------------------------------------------------" << endl; 	
	for(i=0 ; i < nbrOfBasisFunctions ; i++) {		
		phi_a = map.getBasisFunction4Global(i);
		
		for(j=0 ; j < nbrOfBasisFunctions ; j++) {
			
			if(j<i) {
				outFileStream << setw(width) <<"x";
				
			} else {
				phi_b = map.getBasisFunction4Global(j);

				value = overlap.getKineticValue(phi_a,phi_b) / 2.;
				outFileStream << setw(width) << setprecision(precision) << value;
			}	
		}
		outFileStream << endl;
	}
	outFileStream << "---------------------------------------------------------------" << endl; 	
	outFileStream << "-                        End Kinetic                          -" << endl;
	outFileStream << "---------------------------------------------------------------" << endl; 		
	time(&end);
	outFileStream << "- Computation time : " << difftime(end,start) << endl;
	outFileStream << "---------------------------------------------------------------" << endl << endl << endl; 	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function : humanPrintPotential.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void humanPrintPotential(const molecule & mol , ostream & outFileStream , int width=15 , int precision=5) 
{
	int nbrOfBasisFunctions , i, j;
	moleculeMapper map;
	potentialIntegral potential;
	gaussianBasisFunction phi_a , phi_b;
	double value;
	time_t start,end;

	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();

	time(&start);
	outFileStream << "---------------------------------------------------------------" << endl; 	
	outFileStream << "-                      Begin Potential                        -" << endl;
	outFileStream << "---------------------------------------------------------------" << endl; 	
	for(i=0 ; i < nbrOfBasisFunctions ; i++) {		
		phi_a = map.getBasisFunction4Global(i);
		
		for(j=0 ; j < nbrOfBasisFunctions ; j++) {
			
			if(j<i) {
				outFileStream << setw(width) <<"x";
				
			} else {
				phi_b = map.getBasisFunction4Global(j);

				value = 0;
				
				for(int atom = 0 ; atom < mol.getNbrOfAtoms() ; atom++) {
					value -= mol.getNbrOfProtons4Atom(atom) * potential.getPotentialValue(phi_a,phi_b,mol.getPosition4Atom(atom), mol.getDistanceUnit4Atom(atom));					 
				}

				if(fabs(value) <10E-15)
					value = 0;
				outFileStream << setw(width) << setprecision(precision) << value;
			}	
		}
		outFileStream << endl;
	}
	outFileStream << "---------------------------------------------------------------" << endl; 	
	outFileStream << "-                        End Potential                        -" << endl;
	outFileStream << "---------------------------------------------------------------" << endl; 	
	time(&end);
	outFileStream << "- Computation time : " << difftime(end,start) << endl;
	outFileStream << "---------------------------------------------------------------" << endl << endl << endl; 	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function : humanPrintHamiltonian.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void humanPrintHamiltonian(const molecule & mol , ostream & outFileStream , int width = 15 , int precision=5) 
{
	int nbrOfBasisFunctions , i, j;
	moleculeMapper map;
	potentialIntegral potential;
	monoElectronicIntegral overlap;
	gaussianBasisFunction phi_a , phi_b;
	double value;
	time_t start, end;

	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();

	time(&start);
	outFileStream << "---------------------------------------------------------------" << endl; 	
	outFileStream << "-                    Begin Hamiltonian                        -" << endl;
	outFileStream << "---------------------------------------------------------------" << endl; 	
	for(i=0 ; i < nbrOfBasisFunctions ; i++) {		
		phi_a = map.getBasisFunction4Global(i);
		
		for(j=0 ; j < nbrOfBasisFunctions ; j++) {
			
			if(j<i) {
				outFileStream << setw(width) <<"x";
				
			} else {
				phi_b = map.getBasisFunction4Global(j);

				value = overlap.getKineticValue(phi_a,phi_b) / 2.;
;
				
				for(int atom = 0 ; atom < mol.getNbrOfAtoms() ; atom++) {
					value -= mol.getNbrOfProtons4Atom(atom) * potential.getPotentialValue(phi_a,phi_b,mol.getPosition4Atom(atom) , mol.getDistanceUnit4Atom(atom));					 
				}
			
				if(fabs(value) <10E-15)
					value = 0;
				outFileStream << setw(width) <<  setprecision(precision) << value;
			}	
		}
		outFileStream << endl;
	}

	outFileStream << "---------------------------------------------------------------" << endl; 	
	time(&end);
	outFileStream << "- Computation time : " << difftime(end,start) << endl;
	outFileStream << "---------------------------------------------------------------" << endl << endl << endl; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function : humanPrintCoulomb.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void humanPrintCoulomb(const molecule & mol , ostream & outFileStream , int width = 15 , int precision=5) 
{
	int nbrOfBasisFunctions , i, j , k , l;
	moleculeMapper map;
	potentialIntegral potential;
	coulombIntegral coulomb;
	gaussianBasisFunction phi_a , phi_b , phi_c , phi_d;
	double value=0;
	time_t start, end;
	
	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();
	
	time(&start);
	outFileStream << "---------------------------------------------------------------" << endl; 	
	outFileStream << "-                    Begin Coulomb                            -" << endl;
	outFileStream << "---------------------------------------------------------------" << endl; 	

	
		for(i=0 ; i < nbrOfBasisFunctions ; i++) {
			phi_a = map.getBasisFunction4Global(i);				
			
			for(j=0 ; j < nbrOfBasisFunctions ; j++) {
				outFileStream << "Line " << setw(5) << i << " x" << setw(4) << j << " : " << setw(4) << i*nbrOfBasisFunctions + j << " : ";
	
				clog << "Coulomb completed at " << setw(4) << setprecision(2) << 100 * double(i*nbrOfBasisFunctions + j ) / double(nbrOfBasisFunctions * nbrOfBasisFunctions) << endl;
				if(j<i) {
					for(k=0 ; k < nbrOfBasisFunctions ; k++) {		
						for(l=0 ; l < nbrOfBasisFunctions ; l++) {		
							outFileStream << setw(width) <<"x";
						}
					} // fin du for k for l.
					outFileStream << endl;
					continue;
				}	// fin diu if i < j.

				phi_b = map.getBasisFunction4Global(j);
				

				for(k=0 ; k < nbrOfBasisFunctions ; k++) {
				
					if( k < i) {
						for(l=0 ; l < nbrOfBasisFunctions ; l++) {		
							outFileStream << setw(width) <<"x";
						}
						continue;
					} // fin du if k < i.
				
					phi_c = map.getBasisFunction4Global(k);

					for(l=0 ; l < nbrOfBasisFunctions ; l++) {
					
						if( l < k ) {
							outFileStream << setw(width) <<"x";
							continue;
						}

						if(k==i && l<j) {
							outFileStream << setw(width) <<"x";
							continue;
						} 
					
						phi_d = map.getBasisFunction4Global(l);
						value = coulomb.getCoulombValue(phi_a, phi_b , phi_c,phi_d);
						
						if(fabs(value) < 10E-15) {
							value = 0;
						}
						
						outFileStream << setw(width) <<  setprecision(precision) << value;
						
					} // fin du for l.
				} // fin du for k.	
				
				
				outFileStream << endl;
			} // fin du for j.
		} // fin du for i.
				
				
					
	outFileStream << "---------------------------------------------------------------" << endl; 	
	time(&end);
	outFileStream << "- Computation time : " << difftime(end,start) << endl;
	outFileStream << "---------------------------------------------------------------" << endl << endl << endl; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method that prints the matrixes for a given molecule in a human readable way.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void humanPrintMatrixes4Molecule(
																	ostream & outStream , 
																	const molecule & mol , 
																	bool overlap , 
																	bool kinetic , 
																	bool potential , 
																	bool hamilton,
																	bool coulomb)
{

	outStream << "---------------------------------------------------------------" << endl; 	
	outStream << "-                      Begin Molecule                         -" << endl;
	outStream << "---------------------------------------------------------------" << endl; 	
	mol.write(outStream);
	outStream << "---------------------------------------------------------------" << endl; 	
	outStream << "-                      End Molecule                           -" << endl;
	outStream << "---------------------------------------------------------------" << endl  << endl << endl; 	
	
	
	if(overlap) {
		humanPrintOverlap(mol,outStream , 13 , 6);
	} 
	
	if(kinetic) {
		humanPrintKinetic(mol , outStream , 13 , 6);
	} 
	
	if(potential) {
		humanPrintPotential(mol, outStream, 13 , 6);
	} 
	
	if(hamilton) {
		humanPrintHamiltonian(mol, outStream, 13 , 6);
	} 
	
	if(coulomb) {
		humanPrintCoulomb(mol, outStream, 13 , 6);
	}
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function : txtPrintOverlap.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void txtPrintOverlap( const molecule & mol , ostream & outFileStream , int width = 15 , int precision = 5) 
{
	int nbrOfBasisFunctions , i, j;
	moleculeMapper map;
	monoElectronicIntegral overlap;
	gaussianBasisFunction phi_a , phi_b;
	double value;
	time_t start,end;
	
	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();
	
	time(&start);
	
	outFileStream << "# BEGIN OVERLAP" << endl;
	outFileStream << nbrOfBasisFunctions << endl << endl;
	
	for(i=0 ; i < nbrOfBasisFunctions ; i++) {		
		phi_a = map.getBasisFunction4Global(i);
		
		for(j=i ; j < nbrOfBasisFunctions ; j++) {
				
			phi_b = map.getBasisFunction4Global(j);
				
			value = overlap.getOverlapValue(phi_a,phi_b);

			outFileStream << setw(width) << setprecision(precision) << value;
			outFileStream << endl;
		}	
	}
	
	time(&end);
	outFileStream << "# END OVERLAP -  Computation time : " << difftime(end,start) << endl << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function : txtPrintHamiltonian.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void txtPrintKinetic( const molecule & mol , ostream & outFileStream , int width = 15 , int precision = 5) 
{
	
	int nbrOfBasisFunctions , i, j;
	moleculeMapper map;
	monoElectronicIntegral overlap;
	gaussianBasisFunction phi_a , phi_b;
	double value;
	time_t start,end;
	
	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();
	
	time(&start);
	
	outFileStream << "# BEGIN KINETIC" << endl;
	outFileStream << nbrOfBasisFunctions << endl << endl;
	
	for(i=0 ; i < nbrOfBasisFunctions ; i++) {		
		phi_a = map.getBasisFunction4Global(i);
		
		for(j=i ; j < nbrOfBasisFunctions ; j++) {
			
			phi_b = map.getBasisFunction4Global(j);
				
			value = overlap.getKineticValue(phi_a,phi_b) / 2.;
				
			outFileStream << setw(width) << setprecision(precision) << value;
			outFileStream << endl;
		}	
	}
	
	time(&end);
	outFileStream << "# END KINETIC - Computation time : " << difftime(end,start) << endl << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function : txtPrintPotential.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void txtPrintPotential(const molecule & mol , ostream & outFileStream , int width=15 , int precision=5) 
{
	int nbrOfBasisFunctions , i, j;
	moleculeMapper map;
	potentialIntegral potential;
	gaussianBasisFunction phi_a , phi_b;
	double value;
	time_t start,end;
	
	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();
	
	time(&start);
	outFileStream << "# BEGIN POTENTIAL" << endl;
	outFileStream << nbrOfBasisFunctions << endl << endl;
	for(i=0 ; i < nbrOfBasisFunctions ; i++) {		
		phi_a = map.getBasisFunction4Global(i);
		
		for(j=i ; j < nbrOfBasisFunctions ; j++) {
			
				phi_b = map.getBasisFunction4Global(j);
			
				value = 0;				
				for(int atom = 0 ; atom < mol.getNbrOfAtoms() ; atom++) {
					value -= mol.getNbrOfProtons4Atom(atom) * potential.getPotentialValue(phi_a,phi_b,mol.getPosition4Atom(atom) , mol.getDistanceUnit4Atom(atom));					 
				}
				
				outFileStream << setw(width) << setprecision(precision) << value;
				outFileStream << endl;
		}
	}

	time(&end);
	outFileStream << "# END POTENTIAL - Computation time : " << difftime(end,start) << endl << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function : txtPrintHamiltonian.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void txtPrintHamiltonian(const molecule & mol , ostream & outFileStream , int width = 15 , int precision=5) 
{
	int nbrOfBasisFunctions , i, j;
	moleculeMapper map;
	potentialIntegral potential;
	monoElectronicIntegral overlap;
	gaussianBasisFunction phi_a , phi_b;
	double value;
	time_t start, end;
	
	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();
	
	time(&start);
	outFileStream << "# BEGIN HAMILTON" << endl; 	
	
	outFileStream << nbrOfBasisFunctions << endl << endl; 	
	
	for(i=0 ; i < nbrOfBasisFunctions ; i++) {		
		phi_a = map.getBasisFunction4Global(i);
		
		for(j=i ; j < nbrOfBasisFunctions ; j++) {
								
				phi_b = map.getBasisFunction4Global(j);
				
				value = overlap.getKineticValue(phi_a,phi_b) / 2.;
				
				for(int atom = 0 ; atom < mol.getNbrOfAtoms() ; atom++) {
					value -= mol.getNbrOfProtons4Atom(atom) * potential.getPotentialValue(phi_a,phi_b,mol.getPosition4Atom(atom) , mol.getDistanceUnit4Atom(atom));					 
				}
				
				if(fabs(value) <10E-15)
					value = 0;
				outFileStream << setw(width) <<  setprecision(precision) << value;
				outFileStream << endl;
			}	
	}
	
	time(&end);
	outFileStream << "# END HAMILTON - Computation time : " << difftime(end,start) << endl << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function : txtPrintCoulomb.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void txtPrintCoulomb(const molecule & mol , ostream & outFileStream , int width = 15 , int precision=5) 
{
	int nbrOfBasisFunctions , i, j , k , l;
	moleculeMapper map;
	potentialIntegral potential;
	coulombIntegral coulomb;
	gaussianBasisFunction phi_a , phi_b , phi_c , phi_d;
	double value=0;
	time_t start, end;
	
	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();
	
	time(&start);
	outFileStream << "# BEGIN COULOMB"<< endl ;
	outFileStream << nbrOfBasisFunctions << endl << endl;
	
		for(i=0 ; i < nbrOfBasisFunctions ; i++) {
			phi_a = map.getBasisFunction4Global(i);				
			
			for(j=0 ; j < nbrOfBasisFunctions ; j++) {
				
				clog << "Coulomb completed at " << setw(4) << setprecision(2) << 100 * double(i*nbrOfBasisFunctions + j ) / double(nbrOfBasisFunctions * nbrOfBasisFunctions) << endl;

				if(j<i) {
					continue;
				}	// fin du if i < j.
				
				phi_b = map.getBasisFunction4Global(j);
				
				
				for(k=0 ; k < nbrOfBasisFunctions ; k++) {
					
					if( k < i) {
						continue;
					} // fin du if k < i.
					
					phi_c = map.getBasisFunction4Global(k);
					
					for(l=0 ; l < nbrOfBasisFunctions ; l++) {
						
						if( l < k ) {
							continue;
						}
						
						if(k==i && l<j) {
							continue;
						} 
						
						phi_d = map.getBasisFunction4Global(l);
						value = coulomb.getCoulombValue(phi_a, phi_b , phi_c,phi_d);
						
						if(fabs(value) < 10E-15) {
							value = 0;
						}
						
						outFileStream << setw(width) <<  setprecision(precision) << value;
						outFileStream << endl;
						
					} // fin du for l.
				} // fin du for k.	
			} // fin du for j.
		} // fin du for i.
				
		time(&end);
		outFileStream << "# END COULOMB - Computation time : " << difftime(end,start) << endl << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method that prints the matrixes for a given molecule in a human readable way.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void txtPrintMatrixes4Molecule(ostream & outStream , 
														   const molecule & mol , 
																 bool overlap , 
																 bool kinetic , 
																 bool potential , 
																 bool hamilton,
																 bool coulomb)
{	
	if(overlap) {
		txtPrintOverlap(mol,outStream , 13 , 6);
	} 
	
	if(kinetic) {
		txtPrintKinetic(mol , outStream , 13 , 6);
	} 
	
	if(potential) {
		txtPrintPotential(mol, outStream, 13 , 6);
	} 
	
	if(hamilton) {
		txtPrintHamiltonian(mol, outStream, 13 , 6);
	} 
	
	if(coulomb) {
		txtPrintCoulomb(mol, outStream, 13 , 6);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fonction qui écrit la matrices d'overlap dans un fichier XML.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void xmlPrintOverlap( const molecule & mol , ostream & outFileStream) 
{
	int nbrOfBasisFunctions , i, j;
	moleculeMapper map;
	monoElectronicIntegral overlap;
	gaussianBasisFunction phi_a , phi_b;
	double value;
	time_t start,end;

	map.attach(mol);
	nbrOfBasisFunctions = map.getTotalNbrOfBasisFunctions();

	time(&start);
	outFileStream << "<Overlap>" << endl;
	outFileStream << "<Matrix>" << endl;


	for(i=0 ; i < nbrOfBasisFunctions ; i++) {		
		phi_a = map.getBasisFunction4Global(i);
		
		for(j=i ; j < nbrOfBasisFunctions ; j++) {
			
				phi_b = map.getBasisFunction4Global(j);

				value = overlap.getOverlapValue(phi_a,phi_b);
				outFileStream << "\t<Coefficient>" << endl;
				outFileStream << "\t\t<Row>" << i << "</Row>" << endl;
				outFileStream << "\t\t<Column>" << j << "</Column>" << endl;
				outFileStream << "\t\t<Value>" << scientific << value << "</Value>" << endl;
				outFileStream << "\t</Coefficient>" << endl;

		}
	}
	outFileStream << "</Matrix>" << endl;
	time(&end);
	outFileStream << "<ComputationTime unit=\"s\">"<< difftime(end,start) << "</ComputationTime>" << endl;
	outFileStream << "</Overlap>" << endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Method that prints the matrixes for a given molecule in xml.
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void xmlPrintMatrixes4Molecule(
																	ostream & outStream , 
																	const molecule & mol , 
																	bool overlap , 
																	bool kinetic , 
																	bool potential , 
																	bool hamilton,
																	bool coulomb)
{
	
	outStream << "<? xml ?>" << endl;
	outStream << "<ChemicalMatrixes>" << endl;
	
	mol.writeXML(outStream);
	
	if(overlap) {
		xmlPrintOverlap(mol,outStream);
	} 
	
	if(kinetic) {
		//xmlPrintKinetic(mol , outStream , 13 , 6);
	} 
	
	if(potential) {
		//xmlPrintPotential(mol, outStream, 13 , 6);
	} 
	
	if(hamilton) {
		//xmlPrintHamiltonian(mol, outStream, 13 , 6);
	} 
	
	if(coulomb) {
		//xmlPrintCoulomb(mol, outStream, 13 , 6);
	}
	
	outStream << "</ChemicalMatrixes>" << endl; 	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THE MAIN ...
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc , char * argv[]) 
{
	bool coulomb , kinetic , hamilton , overlap , potential;
	
	cmdLineParser cmdArgs(argc,argv);
	molecule mol;
	moleculeDocumentParser moleculeParser;
	ofstream outFileStream;
	string inFileName , outFileName;
	
	if(cmdArgs.hasArgument("help") || cmdArgs.empty()) {
		cerr << "Help : for using computeChemicalMatrix : " << endl;
		cerr << "Help : " << endl;
		cerr << "Help : --coulomb                     : tells that the coulomb tensor should be printed." << endl;
		cerr << "Help : --format < human | xml | txt> : tells that the out put format that is wanted (default is human)." << endl;
		cerr << "Help : --hamilton                    : tells you want the printing of the hamilton matrix." << endl;
		cerr << "Help : --input <moleculeFile.xml>    : the file containing data about the molecule." << endl;
		cerr << "Help : --kinetic                     : tells you want the printing of the kinetic matrix." << endl;
		cerr << "Help : --output <outputFile.txt>     : the file in which the matrix will be written (default is std output)." << endl;
		cerr << "Help : --overlap                     : tells you want the printing of the overlap matrix." << endl;
		cerr << "Help : --potential                   : tells you want the printing of the potential matrix." << endl;
		return 0;
	}
	
	
	if(cmdArgs.hasArgument("input") && !(inFileName = cmdArgs.getArgumentOption("input")).empty()) {
		mol = moleculeParser.getMolecule4Document(inFileName);
		//mol.changeDistanceUnit4Atoms(atom::ATOMIC_UNIT);
	} else {
		cerr << "Error : in computeChemicalMatrix." << endl;
		cerr << "Error : no input file was found." << endl;
		cerr << "Error : computeChemicalMatrix --help for more informations" << endl;
		return 1;
	}
		
	if(cmdArgs.hasArgument("output") && !(outFileName = cmdArgs.getArgumentOption("output")).empty()) {
		outFileStream.open(outFileName.c_str());
	
		if(!outFileStream.is_open()) {
			cerr << "Error : in computeChemicalMatrix" << endl;
			cerr << "Error : unable to open file with name \"" << outFileName << "\"" << endl;
			cerr << "Error : aborting now." << endl;
			return 1;
		}		
	} 
	
	if(cmdArgs.hasArgument("coulomb")) {
		coulomb = true;
	} else {
		coulomb = false;
	}
	
	if(cmdArgs.hasArgument("hamilton")) {
		hamilton = true;
	} else {
		hamilton = false;
	}
	
	if(cmdArgs.hasArgument("kinetic")) {
		kinetic = true;
	} else {
		kinetic = false;
	}

	if(cmdArgs.hasArgument("overlap")) {
		overlap = true;
	} else {
		overlap = false;
	}
	
	if(cmdArgs.hasArgument("potential")) {
	 potential = true;
	} else {
	 potential = false;
	} 
	
	if(cmdArgs.hasArgument("format") && cmdArgs.getArgumentOption("format") == "xml") {
			
		if(outFileStream.is_open()) {
			xmlPrintMatrixes4Molecule(outFileStream , mol , overlap , kinetic , potential , hamilton , coulomb);
		}  else {
			xmlPrintMatrixes4Molecule(cout, mol , overlap , kinetic , potential , hamilton , coulomb);
		}
			
	} else 	if(cmdArgs.hasArgument("format") && cmdArgs.getArgumentOption("format") == "txt") {
	
		if(outFileStream.is_open()) {
			txtPrintMatrixes4Molecule(outFileStream , mol , overlap , kinetic , potential , hamilton , coulomb);
		}  else {
			txtPrintMatrixes4Molecule(cout, mol , overlap , kinetic , potential , hamilton , coulomb);
		}
		
	} else {
		
		if(outFileStream.is_open()) {
			humanPrintMatrixes4Molecule(outFileStream , mol, overlap , kinetic , potential , hamilton , coulomb);
		} else {
			humanPrintMatrixes4Molecule(cout, mol, overlap , kinetic , potential , hamilton , coulomb);
		}
		
	}
	
	
	outFileStream.close();
	return 0;
}

