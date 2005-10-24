//____________________________________________________________________________
/*!

\class    genie::RSListRESPXSec

\brief    Computes the total differential RES cross section RES for a list
          of baryon resonances.

          The cross section is the double differential d^2 xsec/ dQ^2 dW \n
          where \n
          \li Q^2 : momentum transfer ^ 2
          \li W   : invariant mass of the final state hadronic system

          The cross section is computed for an input list of resonances
          as the sum of the Rein-Seghal single resonance cross sections
          weighted with the values of their Breit-Wigner distributions at
          the given W,Q^2 (The user needs to make sure that he does not
          run the single resonance cross section code with a configuration
          that inhibits weighting).

          Concrete implementation of the XSecAlgorithmI interface

\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

*/
//____________________________________________________________________________

#ifndef _RS_RESONANCE_LIST_PARTIAL_XSEC_H_
#define _RS_RESONANCE_LIST_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class RSListRESPXSec : public XSecAlgorithmI {

public:

  RSListRESPXSec();
  RSListRESPXSec(string config);
  virtual ~RSListRESPXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;

private:

  const XSecAlgorithmI * SingleResXSecModel   (void) const;
};

}       // genie namespace

#endif  // _RS_RESONANCE_LIST_PARTIAL_XSEC_H_
