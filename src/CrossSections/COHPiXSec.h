//____________________________________________________________________________
/*!

\class    genie::COHPiXSec

\brief    Computes the cross section for COH neutrino-nucleus pi production.\n
          Is a concrete implementation of the XSecIntegratorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 04, 2004

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COHPI_XSEC_H_
#define _COHPI_XSEC_H_

#include "Base/XSecIntegratorI.h"

namespace genie {

class COHPiXSec : public XSecIntegratorI {
public:
  COHPiXSec();
  COHPiXSec(string config);
  virtual ~COHPiXSec();

  //-- XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //--  members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
};

}       // genie namespace
#endif  // _COHPI_XSEC_H_
