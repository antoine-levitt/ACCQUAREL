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

#include "propagator.h"

// these masses come from NIST website
const double propagator::mass_p = 1836.15267261;
const double propagator::mass_n = 1838.6836598;


///////////////////////////////////////////////////////////////////////////////////
// Default Constructor.
///////////////////////////////////////////////////////////////////////////////////
propagator::propagator(void)
{
  ;
}


///////////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////////
propagator::~propagator(void)
{
  Masses.clear();
  Momenta.clear();
  Forces.clear();
  
}

///////////////////////////////////////////////////////////////////////////////////
// get method for the masses
///////////////////////////////////////////////////////////////////////////////////
const containor<double> & propagator::getMasses(void) const
{
  return Masses;
  
}

///////////////////////////////////////////////////////////////////////////////////
// get method for the momenta
///////////////////////////////////////////////////////////////////////////////////
const matrixFull & propagator::getMomenta(void) const
{
  return Momenta;
  
}

///////////////////////////////////////////////////////////////////////////////////
// get method for the forces
///////////////////////////////////////////////////////////////////////////////////
const matrixFull & propagator::getForces(void) const
{
  return Forces;
  
}

///////////////////////////////////////////////////////////////////////////////////
// Method SET for the masses
///////////////////////////////////////////////////////////////////////////////////
void propagator::setMasses(double mass)
{
  int NbrOfAtoms_ = Masses.getSizes();
  for (int i=0; i<NbrOfAtoms_; i++) {
	Masses[i] = mass;
  }
}


///////////////////////////////////////////////////////////////////////////////////
//   Method SET for the masses
///////////////////////////////////////////////////////////////////////////////////
void propagator::setMasses(containor<double> masses_)
{
  assert(Masses.getSizes() == masses_.getSizes());
  Masses = masses_;
}


///////////////////////////////////////////////////////////////////////////////////
//   Method SET for the momenta
///////////////////////////////////////////////////////////////////////////////////
void propagator::setMomenta(double p)
{
  int NbrOfAtoms_ = Momenta.getNbrOfRows();
  int dim_ = Momenta.getNbrOfColumns();
  for (int i=0; i<NbrOfAtoms_; i++) {
	for (int j=0; j<dim_; j++) {
	  Momenta.setCoefficient(i,j,p);
	}
  }
}


///////////////////////////////////////////////////////////////////////////////////
//   Method SET for the momenta
///////////////////////////////////////////////////////////////////////////////////
void propagator::setMomenta(matrixFull p)
{
  assert(Momenta.getNbrOfRows() == p.getNbrOfRows());
  assert(Momenta.getNbrOfColumns() == p.getNbrOfColumns());

  Momenta = p;
}


///////////////////////////////////////////////////////////////////////////////////
// Method to propagate the system for a time dt
///////////////////////////////////////////////////////////////////////////////////
void propagator::velocityVerlet(restrictedHartreeFockWithDensityWithForces & rhf, optimalDamping4RestrictedHartreeFock & oda, molecule & mol, double dt)
{

  int NbrAtoms_ = mol.getNbrOfAtoms();
  assert(NbrAtoms_ == Masses.getSizes());
  assert(NbrAtoms_ == Forces.getNbrOfRows());
  assert(Forces.getNbrOfColumns() == 3);

  // we will need the null matrix
  matrixSymetric matrixNull;
  matrixNull.setMatrixSize(mol.getNbrOfBasisFunctions());
  matrixNull.setCoefficients(0.);
  
  // compute dt*forces/2
  matrixFull forces_r(NbrAtoms_,3);
  for (int i=0; i<NbrAtoms_; i++) {
	for (int dim=0; dim<3; dim++) {
	  forces_r.setCoefficient(i,dim,0.5*dt*Forces(i,dim));
	}
  }

  // half step on momenta
  Momenta.add(forces_r);

  atom::distanceUnit unit_;
  atom atom_;
  dpoint<3> old_position;
  
  // update positions
  for (int i=0; i<NbrAtoms_; i++) {
	atom_ = mol.getAtom(i);
	unit_ = atom_.getDistanceUnit();
	old_position = atom_.getPosition();
	for (int dim=0;dim<3; dim++) {
	  old_position[dim] += dt*Momenta(i,dim)/Masses[i];
	}
	atom_.setPosition(old_position,unit_);
	mol.setAtom(i,atom_);
  }

  // We now have to recompute the forces
  // we recompute Coulomb tensor and so on with the new nuclei coordinates
  rhf.setMolecularSystem(mol);
  // solve electronic problem
  rhf.setDensity(matrixNull);
  rhf.setDensity(oda.getGroundStateDensity(rhf,false));
  Forces = rhf.getForces(mol);
  
  for (int i=0; i<NbrAtoms_; i++) {
	for (int dim=0; dim<3; dim++) {
	  forces_r.setCoefficient(i,dim,0.5*dt*Forces(i,dim));
	}
  }

  // update momenta
  Momenta.add(forces_r);
  
}

///////////////////////////////////////////////////////////////////////////////////
// Method to initialize properly the masses and set the dimension of the momenta
///////////////////////////////////////////////////////////////////////////////////
void propagator::initialize(const molecule & mol)
{
  int NbrAtoms_ = mol.getNbrOfAtoms();
  Momenta.setMatrixSize(NbrAtoms_,3);
  Forces.setMatrixSize(NbrAtoms_,3);
  Masses.setSizes(NbrAtoms_);

  atom atom_;
  int nb_p, nb_n;
  
  for (int i=0; i<NbrAtoms_; i++) {
	nb_p = mol.getAtom(i).getNbrOfProtons();
	nb_n = mol.getAtom(i).getNbrOfNeutrons();
	// in the reduced units, the electron mass = 1
	Masses[i] = mass_p*nb_p + mass_n*nb_n;
  }
  
}

///////////////////////////////////////////////////////////////////////////////////
// Method to initialize properly the masses, set the momenta at 0 and set the correct dimension for the forces
///////////////////////////////////////////////////////////////////////////////////
void propagator::initialize0(const molecule & mol, const matrixFull & f)
{
  initialize(mol);
  setMomenta(0.);
  assert(Forces.getNbrOfRows() == f.getNbrOfRows());
  assert(Forces.getNbrOfColumns() == f.getNbrOfColumns());
  Forces = f;
}



///////////////////////////////////////////////////////////////////////////////////
// Method to compute the kinetic energy
///////////////////////////////////////////////////////////////////////////////////
double propagator::ComputeKineticEnergy(const molecule & mol)
{

  int NbrAtoms_ = Momenta.getNbrOfRows();
  int dim_ = Momenta.getNbrOfColumns();
  assert (dim_ == 3);
  assert(NbrAtoms_ == mol.getNbrOfAtoms());

  double value = 0.;
  double corr_ = 1./pow(atom::AtomicUnit2Angstroem,2);
  double corr= 1.;
  atom::distanceUnit unit_;
  
  for (int i=0; i<NbrAtoms_; i++) {
	// facteur de correction de l'energie cinetique du au choix de l'unite de distance
	unit_ = mol.getAtom(i).getDistanceUnit();
	if (unit_ == atom::ATOMIC_UNIT) {
	  corr = 1.;
	} else {
	  if (unit_ == atom::ANGSTROEM) {
		corr = corr_;
	  } else {
		assert(0);
		corr = 1.;
	  }
	}
	for (int dim=0; dim<dim_; dim++) {
	  value += 0.5*corr*(pow(Momenta(i,dim),2))/Masses[i];
	}
  }
  return value;
}

  
