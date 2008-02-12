//____________________________________________________________________________
/*!

\class    genie::LHAPDF

\brief    LHAPDF library interface.
          Concrete implementation of the LHAPDF interface.

\author   Anselmo Meregaglia <anselmo.meregaglia@cern.ch>, IPHC Strasbourg
          Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, STFC, Rutherford Lab

\created  January 22, 2008

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LHAPDFLIB_H_
#define _LHAPDFLIB_H_

#include "PDF/PDFModelI.h"

namespace genie {

class LHAPDF : public PDFModelI {

public:

  LHAPDF();
  LHAPDF(string config);
  virtual ~LHAPDF();

  //-- implement PDFModelI interface

  double UpValence   (double x, double q2) const;
  double DownValence (double x, double q2) const;
  double UpSea       (double x, double q2) const;
  double DownSea     (double x, double q2) const;
  double Strange     (double x, double q2) const;
  double Charm       (double x, double q2) const;
  double Bottom      (double x, double q2) const;
  double Top         (double x, double q2) const;
  double Gluon       (double x, double q2) const;
  PDF_t  AllPDFs     (double x, double q2) const;

  //-- override the default "Confugure" implementation 
  //   of the Algorithm interface

  void Configure (const Registry & config);
  void Configure (string config);

private:

  void   SetPDFSetFromConfig (void) const;
};

}         // genie namespace

#endif    // _LHAPDF_SET_MODEL_I_H_
