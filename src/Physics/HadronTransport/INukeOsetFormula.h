/**
 * @brief Formula-based implementation of Oset model
 * 
 * @author Tomasz Golan
 * @date 2015
 * @warning Applicable for pion with Tk < 350 MeV
 * @remarks Based on E. Oset et al., Nucl. Phys. A484 (1988) 557-592
 * 
*/

#ifndef INUKE_OSET_FORMULA_H
#define INUKE_OSET_FORMULA_H

#include "INukeOset.h"

// calculate cross section for piN based on Oset model
class INukeOsetFormula : public INukeOset
{
  public:

  //! use to set up Oset class (assign pion Tk, nuclear density etc)
  void setupOset (const double &density, const double &pionTk,
                  const int &pionPDG, const double &protonFraction);

  private:

  // constants 
  static const double fCouplingConstant; //!< f*^2
  static const double fNucleonMass;      //!< average nucleon mass [MeV]
  static const double fNucleonMass2;     //!< average nucleon mass squared
  static const double fDeltaMass;        //!< delta mass in MeV
  static const double fNormalDensity;    //!< normal nuclear density [fm-3]
  static const double fNormFactor;       //!< MeV^-2 -> mb

  static const double fCoefSigma[fNChannels]; //!< s-wave parametrization eq. 3.7
  static const double fCoefB[fNChannels];     //!< s-wave parametrization eq. 3.8
  static const double fCoefD[fNChannels];     //!< s-wave parametrization eq. 3.8
  static const double ImB0;                 //!< s-wave parametrization eq. 3.12

  static const double    fCoefCQ[fNChannels]; //!< quasi-elastic term (eq. 2.21)
  static const double   fCoefCA2[fNChannels]; //!< two-body absorption (eq. 2.21)
  static const double   fCoefCA3[fNChannels]; //!< three-body absorption (eq. 2.21)
  static const double fCoefAlpha[fNChannels]; //!< alpha (eq. 2.21)
  static const double  fCoefBeta[fNChannels]; //!< beta (eq. 2.21)

  // kinematics variables
  double fPionMass;      //!< pion mass in MeV
  double fPionMass2;     //!< pion mass squared
  double fPionMomentum;  //!< pion momentum in MeV
  double fPionEnergy;    //!< pion total energy in MeV
  double fMomentumCMS;   //!< momentum in CMS in MeV
  double fMomentumCMS2;  //!< momentum in CMS squared
  double fInvariantMass; //!< inv mass = sqrt(s mandelstam) in MeV

  // delta variables
  double fReducedHalfWidth;     //!< reduced delta half width in MeV
  double fSelfEnergyTotal;      //!< total delta self energy in MeV
  double fSelfEnergyAbsorption; //!< abs part of delta self energy in MeV
  double fDeltaPropagator2;     //!< |delta propagator|^2 in MeV-2

  // nucleus variables
  double fFermiMomentum;  //!< Fermi momentum in MeV      
  double fFermiEnergy;    //!< Fermi energy in MeV

  // cross sections
  double fCouplingFactor; //!< (coupling constant / pion mass)^2

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


#endif // INUKE_OSET_FORMULA_H
