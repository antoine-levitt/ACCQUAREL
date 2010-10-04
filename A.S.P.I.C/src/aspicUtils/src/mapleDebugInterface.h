#ifndef _MAPLE_DEBUG_INTERFACE_
#define _MAPLE_DEBUG_INTERFACE_

#include <fstream>
#include <gaussianBasisFunction.h>
#include <hermiteProjector.h>
#include <polynome.h>
#include <iostream>
using namespace std;

///////////////////////////////////////////////////////////////////////////
// Génération d'un fichier maple pour la vérification des calculs.
////////////////////////////////////////////////////////////////////////////
class mapleDebugInterface {

private:
	
	////////////////////////////////////
	// Le flux dans lequel on souhaite
	// écrire les données Mapple.
	////////////////////////////////////
	ofstream OutFileStream;

protected:
	
	///////////////////////////////////////////////////////
	// Méthode pour savoir si le flux sortant est 
	// overt et si l'on peut écrire dedans.
	////////////////////////////////////////////////////////
	bool isOutFileStreamAvailable(void);

public:

	/////////////////////////////////////////////////////
	// Méthode d'accès au flux sortant.
	/////////////////////////////////////////////////////
	ofstream & getOutFileStream(void);

	////////////////////////////////////////////////////////
	// Cette méthode à pour but de construire 
	// des tableaux avec des valeurs dedans dans un fichier
	// mapple.
	/////////////////////////////////////////////////////////
	void initMapleArray(const string & arrayName , const containor<double> & array);

	//////////////////////////////////////////////////////////
	// Cette méthode permet d'écrire dans un fichier Mapple 
	// un ensemblre de procédures qui sont les fonctions de 
	// base gaussiennes.
	///////////////////////////////////////////////////////////
	string initMapleGaussianBasisFunction(const string & functionBaseName , const gaussianBasisFunction & phi);
	
	///////////////////////////////////////////////////////////
	// Méthode qui crée une procédure Maple avec un polynome
	// dedans.
	////////////////////////////////////////////////////////////
	string initMapleGaussianPolynome3D(const string & polyBaseName , const gaussianPolynome3D & poly);

	///////////////////////////////////////////////////////////
	// Méthode qui crée une procédure Maple avec un polynome
	// dedans.
	////////////////////////////////////////////////////////////
	string initMapleHermiteProjector(const string & polyBaseName , const hermiteProjector & poly);

	///////////////////////////////////////////////////////////
	// Méthode qui crée une procédure Maple avec un polynome
	// dedans.
	////////////////////////////////////////////////////////////
	string initMaplePolynome(const string & polyBaseName , const gaussianPolynome3D & poly);

	///////////////////////////////////////////////////////////
	// Méthode qui crée une procédure Maple avec un polynome
	// dedans.
	////////////////////////////////////////////////////////////
	string initMaplePolynome(const string & polyBaseName , const polynome & poly , const double & center =0);

	///////////////////////////////////////////////////////////
	// Méthode qui crée une procédure Maple avec un polynome
	// dedans.
	////////////////////////////////////////////////////////////
	string initMaplePolynome(const string & polyBaseName , const polynome3D & poly , const dpoint<3> & center = dpoint<3>(0));

	/////////////////////////////////////////////////////
	// Méthode qui construit le flux de données sortant
	// vers le fichier dont le path est la valeur de 
	// la variable filePath.
	//
	// Cette méthode s'assure que le flux est bien fermé
	// d'ouvrir un nouveau flux.
	/////////////////////////////////////////////////////
	void initOutFileStream(const string & filePath);
	
	/////////////////////////////////////////////////////
	// Méthode qui construit le flux de données sortant
	// vers le fichier dont le path est la valeur de 
	// la variable filePath.
	//
	// Cette méthode s'assure que le flux est bien fermé
	// d'ouvrir un nouveau flux.
	/////////////////////////////////////////////////////
	void terminateOutFileStream(void);
};

#endif

