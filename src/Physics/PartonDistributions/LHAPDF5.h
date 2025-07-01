//____________________________________________________________________________
/*!

\class    genie::LHAPDF5

\brief    LHAPDF5 library interface.
          Concrete implementation of the PDFModelI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  June 06, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _LHAPDF5_H_
#define _LHAPDF5_H_

#include "Physics/PartonDistributions/PDFModelI.h"
#include <string>
using std::string;

namespace genie {

class LHAPDF5 : public PDFModelI {

public:

  LHAPDF5();
  LHAPDF5(string config);
  virtual ~LHAPDF5();

  // Implement PDFModelI interface

  double UpValence   (double x, double Q2) const;
  double DownValence (double x, double Q2) const;
  double UpSea       (double x, double Q2) const;
  double DownSea     (double x, double Q2) const;
  double Strange     (double x, double Q2) const;
  double Charm       (double x, double Q2) const;
  double Bottom      (double x, double Q2) const;
  double Top         (double x, double Q2) const;
  double Gluon       (double x, double Q2) const;
  PDF_t  AllPDFs     (double x, double Q2) const;

  // Override the default "Confugure" implementation
  // of the Algorithm interface

  void Configure (const Registry & config);
  void Configure (string config);

private:

  void   Initialize          (void) const;
  void   SetPDFSetFromConfig (void) const;
};

}         // genie namespace

#endif    // _LHAPDF5_H_
