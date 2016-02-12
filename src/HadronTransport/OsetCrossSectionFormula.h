/**
 * @brief Formula-based implementation of Oset model
 * 
 * @author Tomasz Golan
 * @date 2015
 * @warning Applicable for pion with Tk < 350 MeV
 * @remarks Based on E. Oset et al., Nucl. Phys. A484 (1988) 557-592
 * 
*/

#ifndef OSET_CROSS_SECTION_FORMULA_H
#define OSET_CROSS_SECTION_FORMULA_H

#include "OsetCrossSection.h"

// calculate cross section for piN based on Oset model
class OsetCrossSectionFormula : public OsetCrossSection
{
  public:

  //! use to set up Oset class (assign pion Tk, nuclear density etc)
  void setupOset (const double &density, const double &pionTk,
                  const int &pionPDG, const double &protonFraction);

  private:

  // constants 
  static const double couplingConstant; //!< f*^2
  static const double nucleonMass;      //!< average nucleon mass [MeV]
  static const double nucleonMass2;     //!< average nucleon mass squared
  static const double deltaMass;        //!< delta mass in MeV
  static const double normalDensity;    //!< normal nuclear density [fm-3]
  static const double normFactor;       //!< MeV^-2 -> mb

  static const double coefSigma[nChannels]; //!< s-wave parametrization eq. 3.7
  static const double coefB[nChannels];     //!< s-wave parametrization eq. 3.8
  static const double coefD[nChannels];     //!< s-wave parametrization eq. 3.8
  static const double ImB0;                 //!< s-wave parametrization eq. 3.12

  static const double    coefCQ[nChannels]; //!< quasi-elastic term (eq. 2.21)
  static const double   coefCA2[nChannels]; //!< two-body absorption (eq. 2.21)
  static const double   coefCA3[nChannels]; //!< three-body absorption (eq. 2.21)
  static const double coefAlpha[nChannels]; //!< alpha (eq. 2.21)
  static const double  coefBeta[nChannels]; //!< beta (eq. 2.21)

  // kinematics variables
  double pionMass;      //!< pion mass in MeV
  double pionMass2;     //!< pion mass squared
  double pionMomentum;  //!< pion momentum in MeV
  double pionEnergy;    //!< pion total energy in MeV
  double momentumCMS;   //!< momentum in CMS in MeV
  double momentumCMS2;  //!< momentum in CMS squared
  double invariantMass; //!< inv mass = sqrt(s mandelstam) in MeV

  // delta variables
  double reducedHalfWidth;     //!< reduced delta half width in MeV
  double selfEnergyTotal;      //!< total delta self energy in MeV
  double selfEnergyAbsorption; //!< abs part of delta self energy in MeV
  double deltaPropagator2;     //!< |delta propagator|^2 in MeV-2

  // nucleus variables
  double fermiMomentum;  //!< Fermi momentum in MeV      
  double fermiEnergy;    //!< Fermi energy in MeV

  // cross sections
  double couplingFactor; //!< (coupling constant / pion mass)^2

  //! set nuclear density and Fermi momentum / energy
  void setNucleus    (const double &density);
  //! do kinematics
  void setKinematics (const double &pionTk, const bool &isPi0);
  //! set up Delta
  void setDelta ();
  
  //! calculate delta self energy 
  void setSelfEnergy ();

  //! calculalte delta width reduction in nuclear medium
  double deltaReduction () const;

  //! calculalte cross sections for each channel
  void setCrossSections ();
};


#endif //OSET_CALCULATIONS_H
