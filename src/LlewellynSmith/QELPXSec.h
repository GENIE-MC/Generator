//____________________________________________________________________________
/*!

\class    genie::QELPXSec

\brief    Computes the differential Quasi Elastic cross section dxsec/dq^2.\n
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      C.H.Llewellyn Smith, Physics Reports (Sect. C of Physics Letters) 3,
          No. 5  (1972) 261-379

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

#ifndef _QEL_PARTIAL_XSEC_H_
#define _QEL_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Base/QELFormFactors.h"

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

  mutable QELFormFactors fFormFactors;
};

}       // genie namespace

#endif  // _QEL_PARTIAL_XSEC_H_
