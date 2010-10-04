#include <atom.h>
#include <coulombIntegral.h>
#include <iostream>
#include <generalizedEigenSolver.h>
#include <molecule.h>
#include <moleculeMapper.h>
#include <roothan.h> 
#include <optimalDampingAlgorithm.h>
#include <restrictedHartreeFockWithDensity.h>
#include <tensorSymetric.h>
#include <xmlBasisDataBaseInterface.h>

using namespace std;

int main(int argc , char ** argv)
{

	double distance = 0.959;
	double angleDeg;
	double angleRad;
	
	molecule water;
	atom hydrogen1, hydrogen2 , oxygen;
	restrictedHartreeFockWithDensity rhfWater;
	matrixSymetric matrixNull;
	
	//////////////////////////////////////////////////////////////
	// Here we construct a water moecule directly from the program.
	////////////////////////////////////////////////////////////////
	
	water.setName("H2O");
	water.setDescription("This is a water molecule to perform the A.S.P.I.C tests.");
	water.setNbrOfAtoms(3);
	
	hydrogen1.setChemicalElement("H");
	water.setBasis4Atom(0,xmlBasisDataBaseInterface::getBasisElement("321","H"));

	oxygen.setChemicalElement("O");
	water.setBasis4Atom(1,xmlBasisDataBaseInterface::getBasisElement("321","O"));

	hydrogen2.setChemicalElement("H");
	water.setBasis4Atom(2,xmlBasisDataBaseInterface::getBasisElement("321","H"));

	matrixNull.setMatrixSize(water.getNbrOfBasisFunctions());
	matrixNull.setCoefficients(0);
		
	for(angleDeg = 108; angleDeg <= 108 ; angleDeg += 2) {
		
		angleRad = (angleDeg * M_PI) / 360.;
		
		hydrogen1.setPosition(distance*cos(angleRad), distance*sin(angleRad),0,atom::ANGSTROEM);
		hydrogen2.setPosition(distance*cos(angleRad),-distance*sin(angleRad),0,atom::ANGSTROEM);
		oxygen.setPosition(0,0,0,atom::ANGSTROEM);
		
		water.setAtom(0,hydrogen1);
		water.setAtom(1,oxygen);
		water.setAtom(2,hydrogen2);
		
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// the we use a RHF model.
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << " -----------------------------------------------------------------------------------------" << endl;
	cout << " Building Restricted Hartree Fock Model. " << endl;
	cout << " -----------------------------------------------------------------------------------------" << endl;
	rhfWater.setMolecularSystem(water);
	rhfWater.setDensity(matrixNull);
	cout << "Computation of the overlap  matrix took: " <<  rhfWater.getOverlapComputationTime() << " second(s)." << endl;
	cout << "Computation of the hamilton matrix took: " <<  rhfWater.getHamiltonComputationTime() << " second(s)." << endl;
	cout << "Computation of the coulomb tensor took: " <<  rhfWater.getCoulombComputationTime() << " second(s)." << endl;
	cout << " -----------------------------------------------------------------------------------------" << endl;
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// We use roothan for optimisation
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	roothan4RestrictedHartreeFock roothan4RHF;
	roothan4RHF.setNbrOfIterationsMax(50);
        rhfWater.setDensity(roothan4RHF.getGroundStateDensity(rhfWater,true));
	cout << " Number of Iterations: " << roothan4RHF.getNbrOfIterations() << endl;
	cout << " Final Energy        : " << rhfWater.getEnergy() << endl;
	cout << " EigenValues         : " << rhfWater.getEigenValues() << endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// We use ODA for optimisation
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	rhfWater.setDensity(matrixNull);
	optimalDamping4RestrictedHartreeFock oda4RHF;
	oda4RHF.setNbrOfIterationsMax(20);	
	rhfWater.setDensity(oda4RHF.getGroundStateDensity(rhfWater,true));
	cout << " Number of Iterations: " << oda4RHF.getNbrOfIterations() << endl;
	cout << " Final Energy        : " << rhfWater.getEnergy() << endl;
	cout << " EigenValues         : " << rhfWater.getEigenValues() << endl;
	}
	
	return 0;
}
