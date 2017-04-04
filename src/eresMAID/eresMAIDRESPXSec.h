//____________________________________________________________________________
/*!

\class    genie::eresMAIDRESPXSec

\brief    Computes the double differential cross section for resonance 
          electro- or neutrino-production according to the Rein-Sehgal model.

          The computed cross section is the d^2 xsec/ dQ^2 dW \n

          where \n
            \li \c Q^2 : momentum transfer ^ 2
            \li \c W   : invariant mass of the final state hadronic system

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      D.Rein and L.M.Sehgal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 05, 2004

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _eresMAID_PXSEC_H_
#define _eresMAID_PXSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResonance.h"
#include "eresMAID/FKR.h"

namespace genie {

class MAIDHelicityAmplModelI;
class Spline;
class XSecIntegratorI;

class eresMAIDRESPXSec : public XSecAlgorithmI {

public:
  eresMAIDRESPXSec();
  eresMAIDRESPXSec(string config);
  virtual ~eresMAIDRESPXSec();

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

  mutable FKR fFKR;

  const MAIDHelicityAmplModelI * fHAmplModelCC;
  const MAIDHelicityAmplModelI * fHAmplModelNCp;
  const MAIDHelicityAmplModelI * fHAmplModelNCn;
  const MAIDHelicityAmplModelI * fHAmplModelEMp;
  const MAIDHelicityAmplModelI * fHAmplModelEMn;

  // configuration data
  bool     fWghtBW;            ///< weight with resonance breit-wigner?
  double   fZeta;              ///< FKR parameter Zeta
  double   fOmega;             ///< FKR parameter Omega
  double   fMa2;               ///< (axial mass)^2
  double   fMv2;               ///< (vector mass)^2
  double   fSin48w;            ///< sin^4(Weingberg angle)
  bool     fUsingDisResJoin;   ///< use a DIS/RES joining scheme?
  bool     fUsingNuTauScaling; ///< use NeuGEN nutau xsec reduction factors?
  double   fWcut;              ///< apply DIS/RES joining scheme < Wcut
  double   fN2ResMaxNWidths;   ///< limits allowed phase space for n=2 res
  double   fN0ResMaxNWidths;   ///< limits allowed phase space for n=0 res
  double   fGnResMaxNWidths;   ///< limits allowed phase space for other res
  Spline * fNuTauRdSpl;        ///< xsec reduction spline for nu_tau
  Spline * fNuTauBarRdSpl;     ///< xsec reduction spline for nu_tau_bar
  double   fXSecScaleCC;       ///< external CC xsec scaling factor
  double   fXSecScaleNC;       ///< external NC xsec scaling factor

  const XSecIntegratorI * fXSecIntegrator;
};

}       // genie namespace

#endif  // _REIN_SEHGAL_RES_PXSEC_H_
