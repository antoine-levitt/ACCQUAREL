#include <cmdLineParser.h>
#include <generalizedEigenSolver.h>
#include <iostream>
#include <matrixSymetric.h>
#include <xmlMatrixDocumentParser.h>
using namespace std;


void printUsage(void) 
{
	cerr << "Error : in generalizedEigenSolverSample.exe" << endl;
	cerr << "Error : argument parsing failed." << endl;
	cerr << "--matrix <matrixFileName.xml>" << endl;
	cerr << "--mass <matrixFileName.xml>" << endl;
	
	exit(1);
}


void printMapleEigenSolver(const matrixSymetric & Matrix , const matrixSymetric & Mass , const generalizedEigenSolver &  Solver, ostream & outStream) 
{
	
	assert(Matrix.getNbrOfRows() == Mass.getNbrOfRows());
	
	int nbrOfRows , i , j;
	
	outStream << ">restart;" << endl;
	outStream << ">with(linalg);"<< endl;
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// First we store in the maple way the A matrix.
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	nbrOfRows = Matrix.getNbrOfRows();
	
	outStream << ">A := matrix(" << nbrOfRows << " , " << nbrOfRows << " , [" << endl;
	
	
	for(i=0 ; i < nbrOfRows ; i++) {
		outStream << ">" ;
		for(j=0 ; j < nbrOfRows ; j++) {
			if(j>0) 
				outStream << " , ";
			outStream << Matrix(i,j);
		}
		
		if( i < (nbrOfRows -1))
			outStream << " ,";
		outStream << endl;
	}
	outStream << "> ]);" << endl;
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// First we store in the maple way the B matrix.
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	nbrOfRows = Mass.getNbrOfRows();
	
	outStream << ">B := matrix(" << nbrOfRows << " , " << nbrOfRows << " , [" << endl;
	
	for(i=0 ; i < nbrOfRows ; i++) {
		outStream << ">" ;
		for(j=0 ; j < nbrOfRows ; j++) {
			if(j>0) 
				outStream << " , ";
			outStream << Mass(i,j);
		}
		
		if( i < (nbrOfRows -1))
			outStream << " ,";
		outStream << endl;
	}
	outStream << "> ]);" << endl;
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Then we store the eigen values found.
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	nbrOfRows = Solver.getNbrOfEigenValues();
	
	outStream << "> EigenValues := array(1 .." << nbrOfRows << " , [ " << endl;
	 
	for(i=0 ; i < nbrOfRows ; i++) {
		if(i>0) {
			outStream << " , ";
		} else {
			outStream << "> ";
		}
		outStream << Solver.getEigenValue(i);
	}

	outStream << endl << "> ]);" << endl;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Then we store the eigen vectors found.
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	for (j=0 ; j < nbrOfRows ; j++) {
	outStream << "> EigenVector[ " << j << " ] := array(1 .." << nbrOfRows << " , [ " << endl;
	for(i=0 ; i < nbrOfRows ; i++) {
		if(i>0) {
			outStream << " , ";
		} else {
			outStream << "> ";
		}
		outStream << Solver.getEigenVector(j)[i];
	}
	
		outStream << endl << "> ]);" << endl;
	}
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Then we store the eigen vectors found as a Matrix. 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	outStream << ">EigenVectors := matrix(" << nbrOfRows << " , " << nbrOfRows << " , [" << endl;
	
	
	for(i=0 ; i < nbrOfRows ; i++) {
		outStream << ">" ;
		for(j=0 ; j < nbrOfRows ; j++) {
			if(j>0) 
				outStream << " , ";
			outStream << Solver.getEigenVector(j)[i];
		}
		
		if( i < (nbrOfRows -1))
			outStream << " ,";
		outStream << endl;
	}
	outStream << "> ]);" << endl;
	
	
}


int main(int argc ,char ** argv) 
{


	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// We use the command line parser to perform the analysis of the cmd line.
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	cmdLineParser cmdLine(argc,argv);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// definition of the objects that are gonna be used in this program.
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	generalizedEigenSolver Solver;
	matrixSymetric Matrix , Mass;
	string matrixFileName;	
	xmlMatrixDocumentParser matrixParser;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// At this point of the file we read the matrix.
	//
	// This matrix is defined in an XML file and given to the program with the --matrix option.
	// if this option does not exists in the cmd line, then we print the usage and quit with signal 1.
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	if( cmdLine.hasArgument("matrix") && !(matrixFileName = cmdLine.getArgumentOption("matrix")).empty())
	{
		Matrix = matrixParser.getMatrixSymetricFull4Document(matrixFileName);
		matrixParser.close();
	}	else {
		printUsage();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// At this point of the file we read the mass matrix.
	//
	// This matrix is defined in an XML file and given to the program with the --mass option.
	// if this option does not exists in the cmd line, then we print the usage and quit with signal 1.
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	if( cmdLine.hasArgument("mass") && !(matrixFileName = cmdLine.getArgumentOption("mass")).empty())
	{
		Mass = matrixParser.getMatrixSymetricFull4Document(matrixFileName);
		matrixParser.close();
	}	else {
		printUsage();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Then we pass the arguments to the generalized eigen solver that will perform the diagonalisation.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	Solver.solve(Matrix , Mass);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Then we create the Maple file that is used to verify the results.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(cmdLine.hasArgument("output") && !(matrixFileName = cmdLine.getArgumentOption("output")).empty()) {
		
		ofstream outFileStream(matrixFileName.c_str()); 
		
		if(!outFileStream.is_open()) {
			cerr << "Error : in generalizedEigenSolverSample" << endl;
			cerr << "Error : unable to open file with name \"" << matrixFileName << "\"" << endl;
			cerr << "Error : aborting now." << endl;
			return 1;
		}		
		
		printMapleEigenSolver(Matrix,Mass,Solver, outFileStream);
													
	} else {
		
		printMapleEigenSolver(Matrix,Mass,Solver,cout);
	
	}
	
}










