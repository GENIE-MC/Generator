//____________________________________________________________________________
/*!

\class    genie::BYPDFModel

\brief    Computes corrected PDFs according to the Bodek-Yang model.

          Concrete implementation of the PDFModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 29, 2004

*/
//____________________________________________________________________________

#ifndef _BODEK_YANG_PDF_H_
#define _BODEK_YANG_PDF_H_

#include "PDF/PDFModelI.h"

namespace genie {

class BYPDFModel : public PDFModelI {

public:

  BYPDFModel();
  BYPDFModel(string config);
  virtual ~BYPDFModel();

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

private:
  double DeltaDU (double x) const;
};

}         // genie namespace

#endif    // _BODEK_YANG_PDF_H_
