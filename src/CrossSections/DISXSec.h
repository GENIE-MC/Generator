//____________________________________________________________________________
/*!

\class    genie::DISXSec

\brief    Computes the DIS Cross Section.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#ifndef _DIS_XSEC_H_
#define _DIS_XSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class IntegratorI;

class DISXSec : public XSecAlgorithmI {

public:

  DISXSec();
  DISXSec(string config);
  virtual ~DISXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfigData           (void);
  void LoadSubAlg               (void);
  bool IsWithinIntegrationRange (const Interaction * interaction) const;

  const XSecAlgorithmI * fPartialXSecAlg;
  const IntegratorI *    fIntegrator;

  int    fNlogx;
  int    fNlogy;
  double fXmin;
  double fXmax;
  double fYmin;
  double fYmax;
  double fLogXmax;
  double fLogXmin;
  double fLogYmax;
  double fLogYmin;
  double fdLogX;
  double fdLogY;
  double fWmin;
  double fWmax;
  double fQ2min;
  double fQ2max;
};

}       // genie namespace

#endif  // _DIS_XSEC_H_
