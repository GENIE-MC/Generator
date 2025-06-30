//____________________________________________________________________________
/*!

\class    genie::Born

\brief    Born level nu-electron cross section.

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC & Harvard University

\created  Dec 8, 2021

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BORN_PXSEC_H_
#define _BORN_PXSEC_H_

#include <TComplex.h>

namespace genie {

class Born {

public:
  Born ();
  virtual ~Born ();

  double GetReAlpha (void) { return falpha.Re(); }
  double PXSecCCR      (double s, double t, double mlin, double mlout);
  double PXSecCCV      (double s, double t, double mlin, double mlout);
  double PXSecCCRNC    (double s, double t, double mlin, double mlout);
  double PXSecCCVNC    (double s, double t, double mlin, double mlout);
  double PXSecNCVnu    (double s, double t, double mlin, double mlout);
  double PXSecNCVnubar (double s, double t, double mlin, double mlout);
  double PXSecPhoton   (double s, double t, double mlout2);
  double PXSecPhoton_T (double s12, double s13, double Q2, double ml2);
  double PXSecPhoton_L (double s12, double s13, double Q2, double ml2);
  double GetS           (double mlin, double Enuin);
  double GetT           (double mlin, double mlout, double s, double costhCM);
  double GetU           (double mlin, double mlout, double s, double t);
  bool   IsInPhaseSpace (double mlin, double mlout, double Enuin, double Enuout);
  double Lambda         (double a, double b, double c);

private:

  double fGw;
  double fGz;

  TComplex falpha;
  TComplex fsw2;
  TComplex fcw2;
  TComplex fmw2c;
  TComplex fmz2c;
  TComplex fgae;
  TComplex fgbe;
  TComplex fgav;

};

}       // genie namespace

#endif  // _BORN_H_
