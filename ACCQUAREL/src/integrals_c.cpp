#include<iostream>
#include<atom.h>
#include<monoElectronicIntegral.h>
#include<potentialIntegral.h>
#include<coulombIntegral.h>
#include"integrals_c.h"

double unrolloverlap(const int nbrofprimitives_a,
                     const double center_a[],
                     const double exponents_a[],
                     const double coefficients_a[],
                     const int monomialdegree_a[],
                     const int nbrofprimitives_b,
                     const double center_b[],
                     const double exponents_b[],
                     const double coefficients_b[],
                     const int monomialdegree_b[]) {
  contractions c;
  gaussianBasisFunction phi_a,phi_b;
  monoElectronicIntegral overlap;
  double value;

  c.setNbrOfContractions(nbrofprimitives_a);
  for (int i=0; i<nbrofprimitives_a; i++) {
     c.setExponent(i,exponents_a[i]); c.setCoefficient(i,coefficients_a[i]);
  }
  phi_a.setContractions(c);
  phi_a.setCenter(dpoint<3>(center_a[0],center_a[1],center_a[2]));
  phi_a.setMonomeDegree(ipoint<3>(monomialdegree_a[0],monomialdegree_a[1],monomialdegree_a[2]));
  c.setNbrOfContractions(nbrofprimitives_b);
  for (int i=0; i<nbrofprimitives_b; i++) {
     c.setExponent(i,exponents_b[i]); c.setCoefficient(i,coefficients_b[i]);
  }
  phi_b.setContractions(c);
  phi_b.setCenter(dpoint<3>(center_b[0],center_b[1],center_b[2]));
  phi_b.setMonomeDegree(ipoint<3>(monomialdegree_b[0],monomialdegree_b[1],monomialdegree_b[2]));
  value=overlap.getOverlapValue(phi_a,phi_b);
  return value;
}

double unrollkinetic(const int nbrofprimitives_a,
                     const double center_a[],
                     const double exponents_a[],
                     const double coefficients_a[],
                     const int monomialdegree_a[],
                     const int nbrofprimitives_b,
                     const double center_b[],
                     const double exponents_b[],
                     const double coefficients_b[],
                     const int monomialdegree_b[]) {
  contractions c;
  gaussianBasisFunction phi_a,phi_b;
  monoElectronicIntegral kinetic;
  double value;

  c.setNbrOfContractions(nbrofprimitives_a);
  for (int i=0; i<nbrofprimitives_a; i++) {
     c.setExponent(i,exponents_a[i]); c.setCoefficient(i,coefficients_a[i]);
  }
  phi_a.setContractions(c);
  phi_a.setCenter(dpoint<3>(center_a[0],center_a[1],center_a[2]));
  phi_a.setMonomeDegree(ipoint<3>(monomialdegree_a[0],monomialdegree_a[1],monomialdegree_a[2]));
  c.setNbrOfContractions(nbrofprimitives_b);
  for (int i=0; i<nbrofprimitives_b; i++) {
     c.setExponent(i,exponents_b[i]); c.setCoefficient(i,coefficients_b[i]);
  }
  phi_b.setContractions(c);
  phi_b.setCenter(dpoint<3>(center_b[0],center_b[1],center_b[2]));
  phi_b.setMonomeDegree(ipoint<3>(monomialdegree_b[0],monomialdegree_b[1],monomialdegree_b[2]));
  value=kinetic.getKineticValue(phi_a,phi_b);
  return value;
}

double unrollderiv(const int nbrofprimitives_a,
                   const double center_a[],
                   const double exponents_a[],
                   const double coefficients_a[],
                   const int monomialdegree_a[],
                   const int nbrofprimitives_b,
                   const double center_b[],
                   const double exponents_b[],
                   const double coefficients_b[],
                   const int monomialdegree_b[],
                   const int dimension) {
  contractions c;
  gaussianBasisFunction phi_a,phi_b;
  monoElectronicIntegral deriv;
  double value;

  c.setNbrOfContractions(nbrofprimitives_a);
  for (int i=0; i<nbrofprimitives_a; i++) {
     c.setExponent(i,exponents_a[i]); c.setCoefficient(i,coefficients_a[i]);
  }
  phi_a.setContractions(c);
  phi_a.setCenter(dpoint<3>(center_a[0],center_a[1],center_a[2]));
  phi_a.setMonomeDegree(ipoint<3>(monomialdegree_a[0],monomialdegree_a[1],monomialdegree_a[2]));
  c.setNbrOfContractions(nbrofprimitives_b);
  for (int i=0; i<nbrofprimitives_b; i++) {
     c.setExponent(i,exponents_b[i]); c.setCoefficient(i,coefficients_b[i]);
  }
  phi_b.setContractions(c);
  phi_b.setCenter(dpoint<3>(center_b[0],center_b[1],center_b[2]));
  phi_b.setMonomeDegree(ipoint<3>(monomialdegree_b[0],monomialdegree_b[1],monomialdegree_b[2]));
  value=deriv.getDerivValue(phi_a,phi_b,dimension);
  return value;
}

double unrollpotential(const int nbrofprimitives_a,
                       const double center_a[],
                       const double exponents_a[],
                       const double coefficients_a[],
                       const int monomialdegree_a[],
                       const int nbrofprimitives_b,
                       const double center_b[],
                       const double exponents_b[],
                       const double coefficients_b[],
                       const int monomialdegree_b[],
                       const double center[]) {
  contractions c;
  gaussianBasisFunction phi_a,phi_b;
  potentialIntegral potential;
  atom::distanceUnit unit=atom::ATOMIC_UNIT;
  double value;

  c.setNbrOfContractions(nbrofprimitives_a);
  for (int i=0; i<nbrofprimitives_a; i++) {
     c.setExponent(i,exponents_a[i]); c.setCoefficient(i,coefficients_a[i]);
  }
  phi_a.setContractions(c);
  phi_a.setCenter(dpoint<3>(center_a[0],center_a[1],center_a[2]));
  phi_a.setMonomeDegree(ipoint<3>(monomialdegree_a[0],monomialdegree_a[1],monomialdegree_a[2]));
  c.setNbrOfContractions(nbrofprimitives_b);
  for (int i=0; i<nbrofprimitives_b; i++) {
     c.setExponent(i,exponents_b[i]); c.setCoefficient(i,coefficients_b[i]);
  }
  phi_b.setContractions(c);
  phi_b.setCenter(dpoint<3>(center_b[0],center_b[1],center_b[2]));
  phi_b.setMonomeDegree(ipoint<3>(monomialdegree_b[0],monomialdegree_b[1],monomialdegree_b[2]));
  value=potential.getPotentialValue(phi_a,phi_b,dpoint<3>(center[0],center[1],center[2]),unit);
  return value;
}

double unrollxderiv(const int nbrofprimitives_a,
                    const double center_a[],
                    const double exponents_a[],
                    const double coefficients_a[],
                    const int monomialdegree_a[],
                    const int nbrofprimitives_b,
                    const double center_b[],
                    const double exponents_b[],
                    const double coefficients_b[],
                    const int monomialdegree_b[],
                    const int dimension1,
                    const int dimension2) {
  contractions c;
  gaussianBasisFunction phi_a,phi_b;
  monoElectronicIntegral xderiv;
  double value;

  c.setNbrOfContractions(nbrofprimitives_a);
  for (int i=0; i<nbrofprimitives_a; i++) {
     c.setExponent(i,exponents_a[i]); c.setCoefficient(i,coefficients_a[i]);
  }
  phi_a.setContractions(c);
  phi_a.setCenter(dpoint<3>(center_a[0],center_a[1],center_a[2]));
  phi_a.setMonomeDegree(ipoint<3>(monomialdegree_a[0],monomialdegree_a[1],monomialdegree_a[2]));
  c.setNbrOfContractions(nbrofprimitives_b);
  for (int i=0; i<nbrofprimitives_b; i++) {
     c.setExponent(i,exponents_b[i]); c.setCoefficient(i,coefficients_b[i]);
  }
  phi_b.setContractions(c);
  phi_b.setCenter(dpoint<3>(center_b[0],center_b[1],center_b[2]));
  phi_b.setMonomeDegree(ipoint<3>(monomialdegree_b[0],monomialdegree_b[1],monomialdegree_b[2]));
  value=xderiv.getXDerivValue(phi_a,phi_b,dimension1,dimension2);
  return value;
}

double unrollcoulomb(const int nbrofprimitives_a,
                     const double center_a[],
                     const double exponents_a[],
                     const double coefficients_a[],
                     const int monomialdegree_a[],
                     const int nbrofprimitives_b,
                     const double center_b[],
                     const double exponents_b[],
                     const double coefficients_b[],
                     const int monomialdegree_b[],
                     const int nbrofprimitives_c,
                     const double center_c[],
                     const double exponents_c[],
                     const double coefficients_c[],
                     const int monomialdegree_c[],
                     const int nbrofprimitives_d,
                     const double center_d[],
                     const double exponents_d[],
                     const double coefficients_d[],
                     const int monomialdegree_d[]) {
  contractions c;
  gaussianBasisFunction phi_a,phi_b,phi_c,phi_d;
  coulombIntegral coulomb;
  double value;

  c.setNbrOfContractions(nbrofprimitives_a);
  for (int i=0; i<nbrofprimitives_a; i++) {
     c.setExponent(i,exponents_a[i]); c.setCoefficient(i,coefficients_a[i]);
  }
  phi_a.setContractions(c);
  phi_a.setCenter(dpoint<3>(center_a[0],center_a[1],center_a[2]));
  phi_a.setMonomeDegree(ipoint<3>(monomialdegree_a[0],monomialdegree_a[1],monomialdegree_a[2]));
  c.setNbrOfContractions(nbrofprimitives_b);
  for (int i=0; i<nbrofprimitives_b; i++) {
     c.setExponent(i,exponents_b[i]); c.setCoefficient(i,coefficients_b[i]);
  }
  phi_b.setContractions(c);
  phi_b.setCenter(dpoint<3>(center_b[0],center_b[1],center_b[2]));
  phi_b.setMonomeDegree(ipoint<3>(monomialdegree_b[0],monomialdegree_b[1],monomialdegree_b[2]));
  c.setNbrOfContractions(nbrofprimitives_c);
  for (int i=0; i<nbrofprimitives_c; i++) {
     c.setExponent(i,exponents_c[i]); c.setCoefficient(i,coefficients_c[i]);
  }
  phi_c.setContractions(c);
  phi_c.setCenter(dpoint<3>(center_c[0],center_c[1],center_c[2]));
  phi_c.setMonomeDegree(ipoint<3>(monomialdegree_c[0],monomialdegree_c[1],monomialdegree_c[2]));
  c.setNbrOfContractions(nbrofprimitives_d);
  for (int i=0; i<nbrofprimitives_d; i++) {
     c.setExponent(i,exponents_d[i]); c.setCoefficient(i,coefficients_d[i]);
  }
  phi_d.setContractions(c);
  phi_d.setCenter(dpoint<3>(center_d[0],center_d[1],center_d[2]));
  phi_d.setMonomeDegree(ipoint<3>(monomialdegree_d[0],monomialdegree_d[1],monomialdegree_d[2]));
  value=coulomb.getCoulombValue(phi_a,phi_b,phi_c,phi_d);
  return value;
}
