//____________________________________________________________________________
/*!

\class    genie::QELPXSec

\brief    Computes the differential Quasi Elastic cross section dxsec/dq^2.\n
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

#ifndef _QEL_PARTIAL_XSEC_H_
#define _QEL_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class QELFormFactorsModelI;

class QELPXSec : public XSecAlgorithmI {

public:
  QELPXSec();
  QELPXSec(string config);
  virtual ~QELPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfigData (void);
  void LoadSubAlg     (void);

  const QELFormFactorsModelI * fFormFactorsModel;
};

}       // genie namespace

#endif  // _QEL_PARTIAL_XSEC_H_
