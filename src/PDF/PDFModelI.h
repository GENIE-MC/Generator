//____________________________________________________________________________
/*!

\class    genie::PDFModelI

\brief    Pure abstract base class. Defines the PDFModelI interface to be
          implemented by any algorithmic class computing Parton Density
          Functions.

          Wrapper classes to existing Parton Density Function
          libraries (PDFLIB, LHAPDF) should also adopt this interface so as
          to be integrated into the GENIE framework.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PDF_MODEL_I_H_
#define _PDF_MODEL_I_H_

#include "Algorithm/Algorithm.h"
#include "PDF/PDFt.h"

namespace genie {

class PDFModelI : public Algorithm {

public:

  virtual ~PDFModelI();

  //-- define PDFModelI interface

  virtual double UpValence   (double x, double q2) const = 0;
  virtual double DownValence (double x, double q2) const = 0;
  virtual double UpSea       (double x, double q2) const = 0;
  virtual double DownSea     (double x, double q2) const = 0;
  virtual double Strange     (double x, double q2) const = 0;
  virtual double Charm       (double x, double q2) const = 0;
  virtual double Bottom      (double x, double q2) const = 0;
  virtual double Top         (double x, double q2) const = 0;
  virtual double Gluon       (double x, double q2) const = 0;
  virtual PDF_t  AllPDFs     (double x, double q2) const = 0;

protected:

  PDFModelI();
  PDFModelI(string name);
  PDFModelI(string name, string config);
};

}         // genie namespace

#endif    // _PDF_MODEL_I_H_
