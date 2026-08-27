//____________________________________________________________________________
/*!

\class    genie::DMRESPXSec

\brief    Computes the double differential cross section for resonance
          DM-production according to the Rein-Sehgal model.

          The computed cross section is the d^2 xsec/ dQ^2 dW \n

          where \n
            \li \c Q^2 : momentum transfer ^ 2
            \li \c W   : invariant mass of the final state hadronic system

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      D.Rein and L.M.Sehgal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

          Changes made for Dark Matter - Zach Orr, Colorado State University

\created  May 05, 2004

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _DM_RES_PXSEC_H_
#define _DM_RES_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Physics/BoostedDarkMatter/XSection/FKRDM.h"

namespace genie {

class RSHelicityAmplModelDMI;
class Spline;
class XSecIntegratorI;

class DMRESPXSec : public XSecAlgorithmI {

public:
  DMRESPXSec();
  DMRESPXSec(string config);
  virtual ~DMRESPXSec();

  // implement the XSecAlgorithmI interface
  double XSec         (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral     (const Interaction * i) const;
  bool   ValidProcess (const Interaction * i) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);

  mutable FKRDM fFKRDM;

  const RSHelicityAmplModelDMI * fHAmplModelDMp;
  const RSHelicityAmplModelDMI * fHAmplModelDMn;

  // configuration data
  double fQchiV;               ///< DM vector charge
  double fQchiA;               ///< DM axial charge
  double fQchiS;               ///< DM scalar charge
  int    fVelMode;             ///< DM (scalar/fermion)
  double fMedMass;             ///< DM mediator mass
  double fgZp;                 ///< DM mediator coupling
  bool     fWghtBW;            ///< weight with resonance breit-wigner?
  bool     fNormBW;            ///< normalize resonance breit-wigner to 1?
  double   fZeta;              ///< FKR parameter Zeta
  double   fOmega;             ///< FKR parameter Omega
  double   fMa2;               ///< (axial mass)^2
  double   fMv2;               ///< (vector mass)^2
  bool     fUsingDisResJoin;   ///< use a DIS/RES joining scheme?
  double   fWcut;              ///< apply DIS/RES joining scheme < Wcut
  double   fN2ResMaxNWidths;   ///< limits allowed phase space for n=2 res
  double   fN0ResMaxNWidths;   ///< limits allowed phase space for n=0 res
  double   fGnResMaxNWidths;   ///< limits allowed phase space for other res
  string fKFTable;             ///< table of Fermi momentum (kF) constants for various nuclei
  bool fUseRFGParametrization; ///< use parametrization for fermi momentum insted of table?
  bool fUsePauliBlocking;      ///< account for Pauli blocking?

    double   fSin48w;            ///< sin^4(Weingberg angle)

//  bool     fUsingNuTauScaling; ///< use NeuGEN nutau xsec reduction factors?
//  Spline * fNuTauRdSpl;        ///< xsec reduction spline for nu_tau
//  Spline * fNuTauBarRdSpl;     ///< xsec reduction spline for nu_tau_bar

  double   fXSecScaleDM;       ///< external DM xsec scaling factor


  const XSecIntegratorI * fXSecIntegrator;
};

}       // genie namespace

#endif  // _DM_RES_PXSEC_H_
