//____________________________________________________________________________
/*!

\class    genie::P33PaschosLalakulichPXSec

\brief    Double differential resonance cross section for P33 according to the 
          Paschos, Lalakulich model.

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      O.Lalakulich and E.A.Paschos, Resonance Production by Neutrinos:
          I. J=3/2 Resonances, hep-ph/0501109

\author   This class is based on code written by the model authors (Olga
          Lalakulich, 17.02.2005). The code was modified to fit into the
          GENIE framework by Costas Andreopoulos.

\created  February 22, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _P33_PASCHOS_LALAKULICH_PARTIAL_XSEC_H_
#define _P33_PASCHOS_LALAKULICH_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class BaryonResDataSetI;

class P33PaschosLalakulichPXSec : public XSecAlgorithmI {

public:

  P33PaschosLalakulichPXSec();
  P33PaschosLalakulichPXSec(string name);
  virtual ~P33PaschosLalakulichPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);

  double Pauli   (double Q2, double W, double MN) const; ///< Pauli suppression for D2
  double Nu      (double Q2, double W, double MN) const; ///< kinematic variables
  double NuStar  (double Q2, double W, double MN) const; ///< ...
  double PPiStar (double W, double MN) const; ///< ...

  const BaryonResDataSetI * fRESDataTable;

  bool   fTurnOnPauliCorrection;
  double fMa;
  double fMv;
  double fCos28c;
};

}       // genie namespace

#endif  // _P33_PASCHOS_LALAKULICH_PARTIAL_XSEC_H_

