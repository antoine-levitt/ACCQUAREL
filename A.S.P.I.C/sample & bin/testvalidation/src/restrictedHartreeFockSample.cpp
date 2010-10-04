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
	molecule test;
	atom element;
	restrictedHartreeFockWithDensity rhftest;
	matrixSymetric matrixNull;

	test.setName("Atom");
	test.setDescription("This is a single atom to perform the A.S.P.I.C tests.");
	test.setNbrOfAtoms(1);

	element.setChemicalElement("Ti");
	test.setBasis4Atom(0,xmlBasisDataBaseInterface::getBasisElement("321","Ti"));
        element.setPosition(0,0,0,atom::ANGSTROEM);

	test.setAtom(0,element);

	matrixNull.setMatrixSize(test.getNbrOfBasisFunctions());
	matrixNull.setCoefficients(0);

	cout << " -----------------------------------------------------------------------------------------" << endl;
	cout << " Building Restricted Hartree Fock Model. " << endl;
	cout << " -----------------------------------------------------------------------------------------" << endl;
	rhftest.setMolecularSystem(test);
	rhftest.setDensity(matrixNull);
	cout << "Computation of the overlap  matrix took : " <<  rhftest.getOverlapComputationTime() << " second(s)." << endl;
	cout << "Computation of the hamilton matrix took : " <<  rhftest.getHamiltonComputationTime() << " second(s)." << endl;
	cout << "Computation of the coulomb tensor took  : " <<  rhftest.getCoulombComputationTime() << " second(s)." << endl;
	cout << " -----------------------------------------------------------------------------------------" << endl;

//	roothan4RestrictedHartreeFock roothan4RHF;
//	roothan4RHF.setNbrOfIterationsMax(40);
//	cout << "Computation of the Ground State Done : " << endl;
//	rhftest.setDensity(roothan4RHF.getGroundStateDensity(rhftest,true));
//	cout << " Final Energy         : " << rhftest.getEnergy() << endl;
//	cout << " EigenValues          : " << rhftest.getEigenValues() << endl;

        rhftest.setDensity(matrixNull);
	optimalDamping4RestrictedHartreeFock oda4RHF;
	oda4RHF.setNbrOfIterationsMax(100);	
        oda4RHF.setMax4Error(.000001);
	cout << " Computation of the Ground State Done : " << endl;
	rhftest.setDensity(oda4RHF.getGroundStateDensity(rhftest,true));
	cout << " Final Energy         : " << rhftest.getEnergy() << endl;
	cout << " EigenValues          : " << rhftest.getEigenValues() << endl;

	return 0;
}
