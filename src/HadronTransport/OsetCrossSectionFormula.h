/*
 * Author: T. Golan
 * Based on: E. Oset et al., Nucl. Phys. A484 (1988) 557-592
 * 
*/

#ifndef OSET_CROSS_SECTION_FORMULA_H
#define OSET_CROSS_SECTION_FORMULA_H

#include "OsetCrossSection.h"

// calculate cross section for piN based on Oset model
class OsetCrossSectionFormula : public OsetCrossSection
{
  public:

  // use to set up Oset class (assign pion Tk, nuclear density etc)
  void setupOset (const double &density, const double &pionTk,
                  const int &pionPDG, const double &protonFraction);

  private:

  // constants 
  static const double couplingConstant; // f*^2
  static const double chargedPionMass;  // charged pion mass in MeV
  static const double neutralPionMass;  // charged pion mass in MeV
  static const double chargedPionMass2; // charged pion mass squared
  static const double neutralPionMass2; // charged pion mass squared
  static const double protonMass;       // proton mass in MeV
  static const double neutronMass;      // neutron mass in MeV
  static const double nucleonMass;      // average nucleon mass
  static const double nucleonMass2;     // average nucleon mass squared
  static const double deltaMass;        // delta mass in MeV
  static const double normalDensity;    // normal nuclear density [fm-3]
  static const double normFactor;       // MeV^-2 -> mb

  // s-wave parametrization
  static const double coefSigma[nChannels]; // see. 3.7
  static const double coefB[nChannels];     // see. 3.8
  static const double coefD[nChannels];     // see. 3.8
  static const double ImB0;                 // see. 3.12

  // delta parametrization coefficients (NuclPhys A468 (1987) 631-652)
  static const double    coefCQ[nChannels]; // quasi-elastic term
  static const double   coefCA2[nChannels]; // two-body absorption
  static const double   coefCA3[nChannels]; // three-body absorption
  static const double coefAlpha[nChannels]; // see 2.21
  static const double  coefBeta[nChannels]; // see 2.21

  // kinematics variables
  const double *pionMass;      // pointer to charged / neutral pion mass
  const double *pionMass2;     // pointer to pion mass ^2
  double pionMomentum;         // pion momentum in MeV
  double pionEnergy;           // pion total energy in MeV
  double momentumCMS;          // momentum in CMS in MeV
  double momentumCMS2;         // momentum in CMS squared
  double invariantMass;        // inv mass = sqrt(s mandelstam) in MeV

  // delta variables
  double reducedHalfWidth;     // reduced delta half width in MeV
  double selfEnergyTotal;      // total delta self energy in MeV
  double selfEnergyAbsorption; // abs part of delta self energy in MeV
  double deltaPropagator2;     // |delta propagator|^2 in MeV-2

  // nucleus variables
  double fermiMomentum;  // in MeV      
  double fermiEnergy;    // in MeV

  // cross sections
  double couplingFactor;       // (coupling constant / pion mass)^2

  // set nuclear density and Fermi momentum / energy
  void setNucleus    (const double &density);
  // do kinematics
  void setKinematics (const double &pionTk, const bool &isPi0);
  // calculate delta width / propagator / self-energy
  void setDelta ();
  
  // calculate delta self energy (eq. 2.21)
  void setSelfEnergy ();

  // delta width reduction in nuclear medium
  double deltaReduction () const;

  // set up cross sections for all channels
  void setCrossSections ();
};


#endif //OSET_CALCULATIONS_H
