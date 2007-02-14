//____________________________________________________________________________
/*!

\class    genie::PDFLIB

\brief    PDFLIB library interface.

          Concrete implementation of the PDFModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 06, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _PDFLIB_H_
#define _PDFLIB_H_

#include "PDF/PDFModelI.h"

namespace genie {

class PDFLIB : public PDFModelI {

public:

  PDFLIB();
  PDFLIB(string config);
  virtual ~PDFLIB();

  //-- impement PDFModelI interface

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

  void   Initialize          (void) const;
  void   SetPDFSetFromConfig (void) const;
};

}         // genie namespace

#endif    // _PDF_SET_MODEL_I_H_
