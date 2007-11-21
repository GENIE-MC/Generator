//____________________________________________________________________________
/*!

\class    genie::ReinSeghalRESPXSec

\brief    Computes the double differential cross section for production of a
          single baryon resonance according to the \b Rein-Seghal model.

          The computed cross section is the d^2 xsec/ dQ^2 dW \n

          where \n
            \li \c Q^2 : momentum transfer ^ 2
            \li \c W   : invariant mass of the final state hadronic system

          If it is specified (at the external XML configuration) the cross
          section can weighted with the value of the resonance's Breit-Wigner
          distribution at the given W. The Breit-Wigner distribution type can
          be externally specified. \n

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 05, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _REIN_SEGHAL_RES_PXSEC_H_
#define _REIN_SEGHAL_RES_PXSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResParams.h"
#include "BaryonResonance/BreitWignerI.h"
#include "ReinSeghal/FKR.h"

namespace genie {

class BreitWignerI;
class BaryonResDataSetI;
class RSHelicityAmplModelI;
class Spline;
class XSecIntegratorI;

class ReinSeghalRESPXSec : public XSecAlgorithmI {

public:
  ReinSeghalRESPXSec();
  ReinSeghalRESPXSec(string config);
  virtual ~ReinSeghalRESPXSec();

  //! XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
  ///////////bool   ValidKinematics (const Interaction * i) const;

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);

  mutable FKR fFKR;
  mutable BaryonResParams fBRP;

  //! configuration data

  const BreitWignerI *         fBreitWigner;
  const BaryonResDataSetI *    fBaryonResDataSet;
  const RSHelicityAmplModelI * fHAmplModelCC;
  const RSHelicityAmplModelI * fHAmplModelNCp;
  const RSHelicityAmplModelI * fHAmplModelNCn;

  bool     fWghtBW;            ///< weight with resonance breit-wigner?
  double   fZeta;              ///< FKR parameter Zeta
  double   fOmega;             ///< FKR parameter Omega
  double   fMa2;               ///< (axial mass)^2
  double   fMv2;               ///< (vector mass)^2
  bool     fUsingDisResJoin;   ///< use a DIS/RES joining scheme?
  bool     fUsingNuTauScaling; ///< use NeuGEN nutau xsec reduction factors?
  double   fWcut;              ///< apply DIS/RES joining scheme < Wcut
  double   fN2ResMaxNWidths;   ///< limits allowed phase space for n=2 res
  double   fN0ResMaxNWidths;   ///< limits allowed phase space for n=0 res
  double   fGnResMaxNWidths;   ///< limits allowed phase space for other res
  Spline * fNuTauRdSpl;        ///< xsec reduction spline for nu_tau
  Spline * fNuTauBarRdSpl;     ///< xsec reduction spline for nu_tau_bar

  const XSecIntegratorI * fXSecIntegrator;
};

}       // genie namespace

#endif  // _REIN_SEGHAL_RES_PXSEC_H_
