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

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NIEVES_QELCC_CROSS_SECTION_H_
#define _NIEVES_QELCC_CROSS_SECTION_H_

#include "Base/XSecAlgorithmI.h"
#include "Base/QELFormFactors.h"
#include "GHEP/GHepRecord.h"
#include "Nuclear/FermiMomentumTable.h"
#include "Nuclear/NuclearModelI.h"

#include <complex>

namespace genie {

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
  double                       fCos8c2;           ///< cos^2(cabbibo angle)

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

  bool   fDoAvgOverNucleonMomentum;    ///< Average cross section over hit nucleon monentum?
  double fEnergyCutOff;                ///< Average only for energies below this cutoff defining 
                                       ///< the region where nuclear modeling details do matter

  // Variables to do integrals 
  std::vector<double> fIntervalFractions;
  std::vector<double> fW;

  //Functions needed to calculate XSec:

  // Calculates values of CN, CT, CL, and imU, and stores them in the provided 
  // variables. If target is not a nucleus, then CN, CN, and CL are all 1.0.
  // r must be in units of fm.
  void CNCTCLimUcalc(TLorentzVector qTildeP4, double M, double r, 
		     bool is_neutrino, bool tgtIsNucleus, int tgt_pdgc,
		     int A, int Z, int N, bool hitNucIsProton,
		     double & CN, double & CT, double & CL,
		     double & imU, double & t0, double & r00) const;

  //Equations to calculate the relativistic Lindhard function for Amunu
  /*std::complex<double> relLindhardIm(double q0gev, double dqgev,
				     double kFngev, double kFpgev,
				     double M, bool isNeutrino) const;*/
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

  // Functions to do integrals
  std::vector<double> integrationSetup(double a,double b,int n) const;
  double integrate(double a,double b,int n,std::vector<double> y) const;
  //  std::complex<double> integrate(double a,double b,int n,
  //			 std::vector<std::complex<double> > y) const;

  //input must be length 4. Returns 1 if input is an even permutation of 0123,
  //-1 if input is an odd permutation of 0123, and 0 if any two elements
  //are equal
  int leviCivita(int input[]) const;

  double LmunuAnumu(const TLorentzVector neutrinoMom,
		    const TLorentzVector inNucleonMom,
		    const TLorentzVector leptonMom,
		    const TLorentzVector outNucleonMom,
		    double M, double r, bool is_neutrino, 
		    bool tgtIsNucleus,
		    int tgt_pdgc, int A, int Z, int N,
		    bool hitNucIsProton) const;

  void SetRunningLepton(GHepRecord * evrec) const;

  // NOTE: THE REMAINING CODE IS FOR TESTING PURPOSES ONLY

  mutable bool                 fPrintData;        ///< print data
  mutable bool                 fCompareNievesTensors;     ///< print tensors
  mutable TString              fTensorsOutFile;   ///< file to print tensors to
  mutable double               fVc,fCoulombFactor,fresult1,fresult2,
    falpha,fZ,frmax,frcurr;
  mutable double               fElocal, fPlocal,fElep,fPlep;
  //mutable double               q2Orig;
  mutable double               rhopStored,rhonStored,rhoStored,rho0Stored;
  //mutable double               fKF1, fKF2;
  /*mutable double               fc0, fPrimeStored;
  mutable double               fVt,fVl,fRelLinReal,fRelLinIm,
                               fRelLinTotReal,fRelLinTotIm;
  mutable double               fl1,fl2,fl3,q2rellin,fkf,fef,fl2im,fl3im;*/
  void CompareNievesTensors(const Interaction* i) const;

};

}       // genie namespace

#endif  
