//____________________________________________________________________________
/*!

\class    genie::P33PaschosLalakulichPXSec

\brief    Double differential resonance cross section d^2xsec / dQ^2 dW
          for P33 according to the Paschos, Lalakuich model.

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      O.Lalakulich and E.A.Paschos, Resonance Production by Neutrinos:
          I. J=3/2 Resonances, hep-ph/0501109

\author   This class is based on code written by the model authors (Olga
          Lalakulich, 17.02.2005). The code was modified to fit into the
          GENIE framework by Costas Andreopoulos.

\created  February 22, 2005

*/
//____________________________________________________________________________

#ifndef _P33_PASCHOS_LALAKULICH_PARTIAL_XSEC_H_
#define _P33_PASCHOS_LALAKULICH_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class P33PaschosLalakulichPXSec : public XSecAlgorithmI {

public:

  P33PaschosLalakulichPXSec();
  P33PaschosLalakulichPXSec(const char * param_set);
  virtual ~P33PaschosLalakulichPXSec();

  //-- XSecAlgorithmI interface implementation

  double XSec (const Interaction * interaction) const;

private:

  double Pauli   (double Q2, double W) const; ///< Pauli suppression for D2

  double Nu      (double Q2, double W) const; ///< kinematic variables
  double NuStar  (double Q2, double W) const; ///< ...
  double PPiStar (double W)            const; ///< ...

  double MA (void) const; ///< gets the registry value or sets default
  double MV (void) const; ///< ...
};

}       // genie namespace

#endif  // _P33_PASCHOS_LALAKULICH_PARTIAL_XSEC_H_

