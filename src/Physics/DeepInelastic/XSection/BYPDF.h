//____________________________________________________________________________
/*!

\class    genie::BYPDF

\brief    Computes corrected PDFs according to the Bodek-Yang model.

          Concrete implementation of the PDFModelI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  September 29, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _BODEK_YANG_PDF_H_
#define _BODEK_YANG_PDF_H_

#include "Physics/PartonDistributions/PDFModelI.h"

namespace genie {

class BYPDF : public PDFModelI {

public:

  BYPDF();
  BYPDF(string config);
  virtual ~BYPDF();

  //! PDFModelI interface implementation
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

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void   LoadConfig (void);
  double DeltaDU    (double x) const;

  //! configuration parameters

  const PDFModelI * fBasePDFModel; ///< base (uncorrected) PDF model

  double fX0;    ///< correction param X0
  double fX1;    ///< correction param X1
  double fX2;    ///< correction param X2
  double fQ2min; ///< min. Q2 for PDF evaluation
};

}         // genie namespace

#endif    // _BODEK_YANG_PDF_H_
