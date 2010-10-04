/* 
 * The propagator class of A.S.P.I.C. 
 * Written by Frederic Legoll nov 06
 *
 * Copyright (C) 2006  Frederic Legoll
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


#ifndef _PROPAGATOR_
#define _PROPAGATOR_

#include <matrixFull.h>
#include <molecule.h>
#include "restrictedHartreeFockWithDensityWithForces.h"
#include "optimalDampingAlgorithm.h"

/**
 * The propagator class enables one to move in time a molecular system, according 
 * to different algorithms.
 */
class propagator
{

 private:

  /**
   * An array containing the masses of the nuclei that compose the molecule.
   */
  containor<double> Masses;

  /**
   * A matrix containing the momenta of the nuclei
   */
  matrixFull Momenta;

    /**
   * A matrix containing the forces on the nuclei
   */
  matrixFull Forces;

 public:

  /**
   * Default constructor
   */
  propagator(void);

  
  /**
   * Destructor
   */
  ~propagator(void);

  /**
   * Proton mass (in electron mass unit)
   */
  static const double mass_p;

  /**
   * Neutron mass (in electron mass unit)
   */
  static const double mass_n;

  /**
   * Method GET for the masses
   */
  const containor<double> & getMasses(void) const;

  /**
   * Method GET for the momenta
   */
  const matrixFull & getMomenta(void) const;

  /**
   * Method GET for the forces
   */
  const matrixFull & getForces(void) const;

  /**
   * Method SET for the masses
   */
  void setMasses(double mass);

    /**
   * Method SET for the masses
   */
  void setMasses(containor<double> masses_);

  /**
   * Method SET for the momenta
   */
  void setMomenta(double p);

    /**
   * Method SET for the momenta
   */
  void setMomenta(matrixFull p);

  /**
   * Method to propagate the system for a time dt according to velocity verlet algorithm (forces are computed according to the rhf model, with the oda electronic optimisation solver) 
   */
  void velocityVerlet(restrictedHartreeFockWithDensityWithForces & rhf, optimalDamping4RestrictedHartreeFock & oda, molecule & mol, double dt);

  /**
   * Method to initialize properly the masses and set the dimension of the momenta and forces matrix
   */
  void initialize(const molecule & mol);

  /**
	 * Method to initialize properly the masses, set the momenta at 0 and set the correct dimension for the forces
   */
  void initialize0(const molecule & mol, const matrixFull & f);

  /**
   * Method to compute the kinetic energy
   */
  double ComputeKineticEnergy(const molecule & mol);
  
  
};

#endif
