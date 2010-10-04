/* 
* The gaussian integrals library of A.S.P.I.C. 
 * Written and directed by François Lodier support.aspic@gmail.com.
 *
 * Copyright (C) 2005  François Lodier
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef _HERMITE_INTEGRATION_BUFFER_
#define _HERMITE_INTEGRATION_BUFFER_

#include <containor3D.h>
#include <gaussianPolynome.h>
#include "hermiteProjector.h"

/**
 * Cette classe calcul des polynome de la forme
 * v^{r / 2} * H_r(sqrt(v)*R*u) [* exp(-v * R^2 * u^2)].
 * 
 * R = le centre.
 * v = l'exposant.
 * u = la variable.
 *
 * Le seul souci est que rien n'inerdit à R d'etre nul.
 * ceci peut donc poser pas mal de soucis.
 */

class hermiteBasisPolynome 
{
private :
	/**
	 * Ensuite on va aussi stocker le centre.
	 */
	double Center;

	/**
	 * Ensuite ib stocke le degré.
	 */
	int Degree;
	
	/**
	 * Ensuit on va stocker l'exposant.
	 */
	double Exponent;

	/**
	 * Bon on va stocker le polynome sous
	 * cette forme la.
	 */
	gaussianPolynome HermiteBasisPolynome;
	
	/**
	 * Pour savoir si l'initialisation à été 
	 * faite correctement.
	 */
	bool HermiteBasisPolynomeIsInitialized;
	/**
	 * Markeur pour savoir si il faut recalculer
	 * le polynome.
	 */
	bool HermiteBasisPolynomeIsUp2Date;


protected :
	
	/**
	 * Méthode pour faire la copie d'un 
	 * objet hermiteBasisPolynome.
	 */
	void copy(const hermiteBasisPolynome & hBP);

	/**
	 * Méthode d'initialisation des polynomes.
	 */
	void initHermiteBasisPolynome(void);

	/**
	 * Méthode pour calculer le polynome.
	 */
	void computeHermiteBasisPolynome(void);

	/**
	 * Méthode spéciale pour calculer le polynome
	 * lorsque l'on choisit un centre nul.
	 */
	void computeHermiteBasisPolynomeNullCenterException(const int & degree);

	/**
	 * Méthode pour calculer le polynome.
	 */
	void computeNextHermiteBasisPolynome(void);

public:

	/**
	 * Constructeur de la classe Hermite Polynome.
	 */
	hermiteBasisPolynome(void);

	/**
	 * Destructeur de la classe Hermite Polynome.
	 */
	~hermiteBasisPolynome(void);

	/**
	 * Méthode pour ajouter un au degré du polynome
	 */
	void addDegree(void);

	/**
	 * Méthode d'accès au centre.
	 */
	const double & getCenter(void) const;

	/**
	 * Méthode pour connaitre le degré du polynome
	 * de hermite stocké dans l'objet.
	 */
	const int & getDegree(void) const;

	/**
	 * Méthode d'accès à l'Exposant.
	 */
	const double & getExponent(void) const;

	/**
	 * Méthode pour accéder au polynome 
	 */
	const polynome & getPolynome(void) const;

	/**
	 * Méthode pour accéder aux coefficient du 
	 * polynome.
	 */
	double getPolynomeCoefficient(const int & degree) const;

	/**
	 * Méthode pour fixer le centre.
	 */
	void setCenter(const double & center);

	/**
	 * Méthode pour fixer le degré du polynome de Hermite
	 * que l'on souhaite connaitre.
	 */
	void setDegree(const int & degree);
	
	/**
	 * Méthode pour fixer l'exposant.
	 */
	void setExponent(const double & exponent);


	/**
	 * Méthode pour écrire le polynome
	 * dans un flux.
	 */
	void write(ostream & out) const;

};

/**
 * Méthode pour écrire la connerie dans un flux.
 */
extern ostream & operator<< (ostream & out , const hermiteBasisPolynome & hBP);



/**
 * Classe qui calcule et stocke  toutes les intégrales de la forme 
 * int_{0}^{1}{(v*u)^{r}H_r(vRu)exp(-v^{2} R^2 u^2) du}
 * avec  H_r(vRu) =   H_{r_y}(v R_x u)   H_{r_y}(v R_y u)   H_{r_z}(v R_z u) 
 * et 0 \leq r_x \leq r_x^max 
 * et 0 \leq r_y \leq r_y^max
 * et 0 \leq r_z \leq r_z^max
 */
class hermiteIntegrationBuffer 
{
private:

	/**
	 * Ceci est la variable R dans l'expression
	 * que nous cherchons à calculer.
	 */
	dpoint<3> Center;

	/**
	 * ---------- DEPRECATED -------------
	 * Variable permettant de mettre
	 * un coefficient constant devant l'intégrale.
	 */
	double Coefficient;

	/**
	 * Ceci est la variable v dans l'expression 
	 * que nous cherchons à calculer.
	 */
	double Exponent;
	
	/**
	 * Tableau à trois indices qui contient les valeurs
	 * des intégrales.
	 */
	containor3D<double> IntegrationBuffer;
	
	/**
	 * booleen pour qui permet de savoir si le tableau de résultats
	 * doit etre recaclulé.
	 */
	bool IntegrationBufferIsUp2Date;

protected:

	/**
	 * Méthode pour mettre le buffer à zéro.
	 */
	void clearMonomeIntegrationBuffer(void);

	/**
	 * Méthode pour calculer les valeurs du Buffer.
	 */
	void computeHermiteIntegrationBuffer(void);

	/**
	 * Exception : lorsque le buffer ne 
	 * peut pas recevoir de données.
	 */
	void nothing2DoException(void);

	/** ---------- DEPRECATED -------------
	 * Méthode qui permet de savoir si nous
	 * sommes dans le cas particulier où les
	 * deux centres sont identiques.
	 */
	//void sameCenterException(void);

	/** ---------- DEPRECATED -------------
	 * Méthode qui calcule le centre.
	 */
	//void setCenter(const hermiteProjector & hermitePolynome_1 , const hermiteProjector & hermitePolynome_2);

	/** ---------- DEPRECATED -------------
	 * Méthode qui calcule le coefficient
	 * pour l'intégration ...
	 */
	//void setCoefficient(const hermiteProjector & hermitePolynome_1 , const hermiteProjector & hermitePolynome_2);

	/** ---------- DEPRECATED -------------
	 * Méthode pour mettre le buffer à la taille.
	 */
	//void setIntegrationBufferSize(const hermiteProjector & hermitePolynome_1 , const hermiteProjector & hermitePolynome_2);

	/** ---------- DEPRECATED -------------
	 * Méthode pour calculer l'exposant.
	 */
	//void setExponent(const hermiteProjector & hermitePolynome_1 , const hermiteProjector & hermitePolynome_2);
	
public:

	/**
	 * Constructeur de la classe.
	 */
	hermiteIntegrationBuffer(void);

	/**
	 * Méthode pour connaitre le degré maximal
	 * auquel on peut accéder.
	 */
	const ipoint<3> getIntegrationBufferSize(void) const;

	/**
	 * Méthode pour connaitre le centre.
	 */
	const dpoint<3> & getCenter(void) const;

	/**
	 * Méthode pour connaitre le degré maximal
	 * auquel on peut accéder.
	 */
	int getIntegrationBufferSize(const int & dimension) const;
	
	/** ------ DEPRECATED ------------------
	 * Méthode pour connaitre le coefficient.
	 */
	//double getCoefficient(void) const;

	/**
	 * Méthode pour connaitre l'exposant.
	 */
	double getExponent(void) const;
	
	/**
	 * Méthode d'accès aux intégrales.
   */
	const double & getIntegrationValue(const int & rx ,const int & ry ,const int & rz) const;
	
	/**
	 * Méthode d'accès aux intégrales.
   */
	const double & getIntegrationValue(const ipoint<3> & r) const;

	/**
	 * Méthode pour mettre le buffer à la taille.
	 */
	void setIntegrationBufferSize(ipoint<3> degreeMax);
	
	/**
	 * Méthode pour fixer l'exposant.
	 */
	void setExponent(double exponent);
	
	/**
		* Méthode pour fixer le centre.
	 */
	void setCenter(const dpoint<3> & center);
	
	
	/** ------ DEPRECATED ------------------
	 * Méthode pour construire le buffer associé aux projections
	 * de hermite.
	 */
	//void set(const hermiteProjector & hermitePolynome_1 , const hermiteProjector & hermitePolynome_2);
};

#endif 


