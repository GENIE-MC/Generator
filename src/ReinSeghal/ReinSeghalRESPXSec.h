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
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

____________________________________________________________________________*/

#ifndef _REIN_SEGHAL_RES_PARTIAL_XSEC_H_
#define _REIN_SEGHAL_RES_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BreitWignerI.h"

namespace genie {

class BreitWignerI;
class BaryonResDataSetI;
class RSHelicityAmplModelI;

class ReinSeghalRESPXSec : public XSecAlgorithmI {

public:

  ReinSeghalRESPXSec();
  ReinSeghalRESPXSec(string config);
  virtual ~ReinSeghalRESPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadSubAlg     (void);
  void LoadConfigData (void);

  bool   fWghtBW;
  double fZeta;
  double fOmega;
  double fMa2;
  double fMv2;

  const BreitWignerI *         fBreitWigner;
  const BaryonResDataSetI *    fBaryonResDataSet;
  const RSHelicityAmplModelI * fHAmplModel;
};

}       // genie namespace

#endif  // _REIN_SEGHAL_RES_PARTIAL_XSEC_H_
