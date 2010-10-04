#include <aspicConfiguration.h>
#include <monoElectronicIntegral.h>
#include <fstream>
#include <xmlBasisDataBaseInterface.h>

using namespace std;






int main (int argc , char ** argv) 
{

	///////////////////////////////////////////////////////////////////
	// The argument is the name of the basis to normalize.
	///////////////////////////////////////////////////////////////////
	if(argc != 2) {
		cerr << "Error : in basisNormalizer." << endl;
		cerr << "Error : usage is basisNormalizer [basisName]" << endl;
		cerr << "Error : aborting now" << endl;
		exit(1);
	}


	monoElectronicIntegral overlap;
	ofstream out;
	basisElement basisElement;
	gaussianBasisFunction basisFunction;
	int atom , contraction ,  nbrOfAtoms , nbrOfContractions , nbrOfShells , shell;
	string normalizedBasisFunctionFileName;

	xmlBasisDataBaseInterface::connect(argv[1]);

	if( ! xmlBasisDataBaseInterface::connected()) {
		cerr << "Error : in basisNormalizer." << endl;
		cerr << "Error : no basis with name " << argv[1] << " was found." << endl;
		cerr << "Error : aborting now" << endl;
		exit(1);
	}
	
	
	normalizedBasisFunctionFileName = aspicConfiguration::getBasisDataBaseRoot();
	normalizedBasisFunctionFileName += "/";
	normalizedBasisFunctionFileName += argv[1];
	normalizedBasisFunctionFileName += "-n.xml";

	out.open(normalizedBasisFunctionFileName.c_str());

	if(!out.is_open()) {
		cerr << "Error : in basisNormalizer." << endl;
		cerr << "Error : unable to open out put file with name \"" << normalizedBasisFunctionFileName << "\"." << endl;
		cerr << "Error : aborting now" << endl;
		exit(1);
	}


	out << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << endl;
	out << "<GaussianBasis xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">" << endl << endl;


	nbrOfAtoms = xmlBasisDataBaseInterface::getNbrOfBasisElements();


	for(atom=0 ; atom < nbrOfAtoms ; atom++) {
		
		basisElement = xmlBasisDataBaseInterface::getBasisElement(atom);
		out << "<!--" << endl;
		out << "\t Basis element " << basisElement.getBasisElementKey() << "." << endl;
		out << "-->" << endl;
		out << "<BasisElement ElementKey=\"" << basisElement.getBasisElementKey() << "\">" << endl;
		
		
		nbrOfShells = basisElement.getNbrOfShells();
		
		out << "\t<ShellsList>" << endl;
		for( shell=0 ; shell < nbrOfShells ; shell++) {
	
			out << "\t\t<Shell>" << endl;
			out << "\t\t\t<ShellTypeKey>" << basisElement.getShellTypeKey(shell) << "</ShellTypeKey>" << endl;
			
			basisFunction.setCenter(0);
			basisFunction.setContractions(basisElement.getContractions4Shell(shell));
			basisFunction.setMonomeDegree(basisElement.getMonomeDegree(shell,0));
				
			overlap.normalize(basisFunction);		
		

			cout << " Overlap value is now  : "  << setprecision(16) <<  overlap.getOverlapValue(basisFunction,basisFunction) << endl;

			nbrOfContractions = basisFunction.getNbrOfContractions();

			out << "\t\t\t<ContractionsList>" << endl;
			for(contraction = 0 ; contraction < nbrOfContractions ; contraction++) {
				out << "\t\t\t\t<Contraction>" << endl;
				out << "\t\t\t\t\t<Coefficient>" << setprecision(16) << basisFunction.getCoefficient(contraction) << "</Coefficient>" << endl;
				out << "\t\t\t\t\t<Exponent>" << setprecision(16) << basisFunction.getExponent(contraction) << "</Exponent>" << endl;
				out << "\t\t\t\t</Contraction>" << endl;
			}
			out << "\t\t\t</ContractionsList>" << endl;

			out << "\t\t</Shell>" << endl;		
		}
		out << "\t</ShellsList>" << endl;
		out << "</BasisElement>" << endl << endl;

	}

	out << "</GaussianBasis>" << endl;
	return 0;
}

