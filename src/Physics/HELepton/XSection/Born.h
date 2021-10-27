//____________________________________________________________________________
/*!

\class    genie::Born

\brief    Born level nu-electron cross section.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
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
  double PXSecNCV      (double s, double t, double mlin, double mlout);
  double PXSecCCRNC    (double s, double t, double mlin, double mlout);
  double PXSecCCVNC    (double s, double t, double mlin, double mlout);
  double PXSecPhoton   (double s, double t, double mlout2);
  double PXSecPhoton_T (double s12, double s13, double Q2, double ml2);
  double PXSecPhoton_L (double s12, double s13, double Q2, double ml2);
  double GetS           (double mlin, double Enuin);
  double GetT3          (double mlin, double mlout, double s, double costhCM);
  double GetT4          (double mlin, double mlout, double s, double costhCM);
  double GetU           (double mlin, double mlout, double s, double t);
  double GetELab3       (double mlin, double mlout, double u );
  double GetELab4       (double mlin, double mlout, double t );
  bool   IsInPhaseSpace (double mlin, double mlout, double Enuin, double Enuout);

private:

  double Lambda(double a, double b, double c);

  double fGw;
  double fGz;

  TComplex falpha;
  TComplex fsw2;
  TComplex fcw2;
  TComplex fmw2c;
  TComplex fmz2c;
  TComplex fgLe;
  TComplex fgRe;
  TComplex fgLnu;

};

}       // genie namespace

#endif  // _BORN_H_
