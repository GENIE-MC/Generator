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

class ReinSeghalRESPXSec : public XSecAlgorithmI {

public:

  ReinSeghalRESPXSec();
  ReinSeghalRESPXSec(const char * param_set);
  virtual ~ReinSeghalRESPXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;

private:

  const BreitWignerI * BreitWignerAlgorithm(Resonance_t rid) const;

};

}       // genie namespace

#endif  // _REIN_SEGHAL_RES_PARTIAL_XSEC_H_
