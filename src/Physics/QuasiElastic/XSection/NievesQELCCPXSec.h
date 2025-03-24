//____________________________________________________________________________
/*!

\class    genie::NievesQELCCPXSec

\brief    Computes neutrino-nucleon(nucleus) QELCC differential cross section
          with RPA corrections
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      1. PRC70(2004)055503
          2. EPJA25(2005)299-318
          3. Phys.Rept.188(1990)79
          4. CJP46(1996)0673-0720

\author   Joe Johnston, University of Pittsburgh
          Steven Dytman, University of Pittsburgh
          Igor Kakorin, JINR

\created  April 2016

\updated  March 2025

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
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
#include "Physics/Common/QvalueShifter.h"

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
  const  TVector3 & FinalLeptonPolarization (const Interaction* i) const;
  double IntegratedOverMomentum (const Interaction* i, double r, int mod) const;

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

  double                       fXSecCCScale;        ///< external xsec scaling factor for CC
  double                       fXSecNCScale;        ///< external xsec scaling factor for NC
  double                       fXSecEMScale;        ///< external xsec scaling factor for EM
  const QvalueShifter *        fQvalueShifter ;   ///< Optional algorithm to retrieve the qvalue shift for a given target

  double                       fhbarc;            ///< hbar*c in GeV*fm

  // mutable for testing purposes only!
  mutable bool                 fRPA;              ///< use RPA corrections
  bool                         fCoulomb;          ///< use Coulomb corrections

  const NuclearModelI*         fNuclModel;        ///< Nuclear Model for integration
  // Detect whether the nuclear model is local Fermi gas, and store
  // the relativistic Fermi momentum table if not
  bool                         fLFG;
  const FermiMomentumTable *   fKFTable;
  string                       fKFTableName;
  string                       fLindhardFunction;

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

  /// Scaling factor for the Coulomb potential
  double fCoulombScale;

  /// Enum variable describing which method of computing Rmax should be used
  /// for integrating the Coulomb potential
  Nieves_Coulomb_Rmax_t fCoulombRmaxMode;

  //Functions needed to calculate XSec:

  // Calculates values of CN, CT, CL, and imU, and stores them in the provided
  // variables. If target is not a nucleus, then CN, CN, and CL are all 1.0.
  // r must be in units of fm.
  void CNCTCLimUcalc(TLorentzVector qTildeP4, double M, double r, 
    bool tgtIsNucleus, int A, int Z, int N, double & CN, double & CT, 
    double & CL, bool assumeFreeNucleon) const;

  //Relativistic Lindhard function as is in Fortran code provided by j.Nieves
  // Ref.1, Eq.B2
  double relLindhardIm(double q0, double dq, double kFn, double kFp, double M, bool isNeutrino) const;
  // Ref.2, Eq.61
  double ruLinRelX(double q0, double qm, double kf, double m) const;
   
  std::complex<double> LindhardNuclear(double q0, double dq, double kF, double M) const;
  std::complex<double> LindhardDelta  (double q0, double dq, double kF, double M, double rho) const;
                     
  // Potential for coulomb correction
  double vcr(const Target * target, double r) const;
  
  double MaximalRadius(const Target * target) const;
  
  inline int g(int a, int b) const///< metric g^{ab}=g_{ab}=diag(1,-1,-1,1)
  {
      return (a==b)*(2*(a==0) - 1);
  }
  inline int e(int a, int b, int c, int d) const ///< Levi-Chevita symbol, where e_{0123}=+1
  {
      return (b - a)*(c - a)*(d - a)*(c - b)*(d - b)*(d - c)/12;
  }

  double LmunuAnumu(const TLorentzVector neutrinoMom,
    const TLorentzVector inNucleonMom, const TLorentzVector leptonMom,
    const TLorentzVector outNucleonMom, double M, bool is_neutrino,
    const Target& target, bool assumeFreeNucleon) const;
    

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
    
    class NievesQELSmithMonizIntegrand : public ROOT::Math::IBaseFunctionOneDim
    {
     public:
      NievesQELSmithMonizIntegrand(const NievesQELCCPXSec* alg_, const Interaction* interaction_, int mod_);
      ~NievesQELSmithMonizIntegrand();
       // ROOT::Math::IBaseFunctionOneDim interface
       unsigned int                      NDim   (void)       const;
       double                            DoEval (double rin) const;
       ROOT::Math::IBaseFunctionOneDim * Clone  (void)       const;
     private:
        const NievesQELCCPXSec* alg;
        const Interaction* interaction;
        int mod;
    };

   } // wrap namespace
  } // gsl namespace
 } // utils namespace
} // genie namespace


#endif
