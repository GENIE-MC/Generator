//____________________________________________________________________________
/*!

\class    genie::PDFModelI

\brief    Pure abstract base class. Defines the PDFModelI interface to be
          implemented by wrapper classes to existing Parton Density Function
          libraries (PDFLIB, LHAPDF), or by built-in implementations.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PDF_MODEL_I_H_
#define _PDF_MODEL_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Physics/PartonDistributions/PDFt.h"

namespace genie {

class PDFModelI : public Algorithm {

public:

  virtual ~PDFModelI();

  //-- define PDFModelI interface

  virtual double UpValence   (double x, double Q2) const = 0;
  virtual double DownValence (double x, double Q2) const = 0;
  virtual double UpSea       (double x, double Q2) const = 0;
  virtual double DownSea     (double x, double Q2) const = 0;
  virtual double Strange     (double x, double Q2) const = 0;
  virtual double Charm       (double x, double Q2) const = 0;
  virtual double Bottom      (double x, double Q2) const = 0;
  virtual double Top         (double x, double Q2) const = 0;
  virtual double Gluon       (double x, double Q2) const = 0;
  virtual PDF_t  AllPDFs     (double x, double Q2) const = 0;

protected:

  PDFModelI();
  PDFModelI(string name);
  PDFModelI(string name, string config);
};

}         // genie namespace

#endif    // _PDF_MODEL_I_H_
