//____________________________________________________________________________
/*!

\class    genie::NievesQELCCPXSec

\brief    Computes neutrino-nucleon(nucleus) QELCC differential cross section
          with RPA corrections
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      Physical Review C 70, 055503 (2004)

\author   Joe Johnston, University of Pittsburgh
          Steven Dytman, University of Pittsburgh

\created  April 2016

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NIEVES_QELCC_CROSS_SECTION_H_
#define _NIEVES_QELCC_CROSS_SECTION_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include <complex>
#include <Math/IFunction.h>
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/PauliBlocker.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"

namespace genie {

typedef enum EQELRmax {
  // Use the same maximum radius as VertexGenerator (3*R0*A^(1/3))
  kMatchVertexGeneratorRmax,

  // Use the method for calculting Rmax from Nieves' Fortran code
  kMatchNieves
} Nieves_Coulomb_Rmax_t;

class QELFormFactorsModelI;
class XSecIntegratorI;

class NievesQELCCPXSec : public XSecAlgorithmI {

public:
  NievesQELCCPXSec();
  NievesQELCCPXSec(string config);
  virtual ~NievesQELCCPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void LoadConfig (void);

  mutable QELFormFactors       fFormFactors;      ///<
  const QELFormFactorsModelI * fFormFactorsModel; ///<
  const XSecIntegratorI *      fXSecIntegrator;   ///<
  double                       fCos8c2;           ///< cos^2(cabibbo angle)

  double                       fXSecScale;        ///< external xsec scaling factor

  double                       fhbarc;            ///< hbar*c in GeV*fm

  // mutable for testing purposes only!
  mutable bool                         fRPA;              ///< use RPA corrections
  bool                         fCoulomb;          ///< use Coulomb corrections

  const NuclearModelI*         fNuclModel;        ///< Nuclear Model for integration
  // Detect whether the nuclear model is local Fermi gas, and store
  // the relativistic Fermi momentum table if not
  bool                         fLFG;
  const FermiMomentumTable *   fKFTable;
  string                       fKFTableName;

  /// Enum specifying the method to use when calculating the binding energy of
  /// the initial hit nucleon during spline generation
  QELEvGen_BindingMode_t fIntegralNucleonBindingMode;

  /// Cutoff lab-frame probe energy above which the effects of Fermi motion and
  /// binding energy are ignored when computing the total cross section
  double fEnergyCutOff;

  /// Whether to apply Pauli blocking in XSec()
  bool fDoPauliBlocking;
  /// The PauliBlocker instance to use to apply that correction
  const genie::PauliBlocker* fPauliBlocker;

  /// Nuclear radius parameter r = R0*A^(1/3) used to compute the
  /// maximum radius for integration of the Coulomb potential
  /// when matching the VertexGenerator method
  double fR0;

  /// Enum variable describing which method of computing Rmax should be used
  /// for integrating the Coulomb potential
  Nieves_Coulomb_Rmax_t fCoulombRmaxMode;

  //Functions needed to calculate XSec:

  // Calculates values of CN, CT, CL, and imU, and stores them in the provided
  // variables. If target is not a nucleus, then CN, CN, and CL are all 1.0.
  // r must be in units of fm.
  void CNCTCLimUcalc(TLorentzVector qTildeP4, double M, double r,
    bool is_neutrino, bool tgtIsNucleus, int tgt_pdgc, int A, int Z, int N,
    bool hitNucIsProton, double & CN, double & CT, double & CL, double & imU,
    double & t0, double & r00, bool assumeFreeNucleon) const;

  //Equations to calculate the relativistic Lindhard function for Amunu
  std::complex<double> relLindhardIm(double q0gev, double dqgev,
				     double kFngev, double kFpgev,
				     double M, bool isNeutrino,
				     double & t0, double & r00) const;
  std::complex<double> relLindhard(double q0gev, double dqgev,
				   double kFgev, double M,
				   bool isNeutrino,
				   std::complex<double> relLindIm) const;
  std::complex<double> ruLinRelX(double q0, double qm,
				 double kf, double m) const;
  std::complex<double> deltaLindhard(double q0gev, double dqgev,
				     double rho, double kFgev) const;

  // Potential for coulomb correction
  double vcr(const Target * target, double r) const;

  //input must be length 4. Returns 1 if input is an even permutation of 0123,
  //-1 if input is an odd permutation of 0123, and 0 if any two elements
  //are equal
  int leviCivita(int input[]) const;

  double LmunuAnumu(const TLorentzVector neutrinoMom,
    const TLorentzVector inNucleonMom, const TLorentzVector leptonMom,
    const TLorentzVector outNucleonMom, double M, bool is_neutrino,
    const Target& target, bool assumeFreeNucleon) const;

  // NOTE: THE FOLLOWING CODE IS FOR TESTING PURPOSES ONLY
  // Used to print tensor elements and various inputs for comparison to Nieves'
  // fortran code
  mutable bool                 fCompareNievesTensors;     ///< print tensors
  mutable TString              fTensorsOutFile;   ///< file to print tensors to
  mutable double               fVc,fCoulombFactor;
  void CompareNievesTensors(const Interaction* i) const;
  // END TESTING CODE
};
}       // genie namespace

//____________________________________________________________________________
/*!
\class    genie::utils::gsl::wrap::NievesQELIntegrand

\brief    Auxiliary scalar function for integration over the nuclear density
          when calculaing the Coulomb correction in the Nieves QEL xsec model

\author   Joe Johnston, University of Pittsburgh
          Steven Dytman, University of Pittsburgh

\created  June 03, 2016
*/
//____________________________________________________________________________

namespace genie {
 namespace utils {
  namespace gsl   {
   namespace wrap   {

    class NievesQELvcrIntegrand : public ROOT::Math::IBaseFunctionOneDim
    {
     public:
      NievesQELvcrIntegrand(double Rcurr, int A, int Z);
      ~NievesQELvcrIntegrand();
       // ROOT::Math::IBaseFunctionOneDim interface
       unsigned int                      NDim   (void)       const;
       double                            DoEval (double rin) const;
       ROOT::Math::IBaseFunctionOneDim * Clone  (void)       const;
     private:
       double fRcurr;
       double fA;
       double fZ;
    };

   } // wrap namespace
  } // gsl namespace
 } // utils namespace
} // genie namespace


#endif
