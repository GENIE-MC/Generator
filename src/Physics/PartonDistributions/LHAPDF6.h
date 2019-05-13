//____________________________________________________________________________
/*!

\class    genie::LHAPDF6

\brief    LHAPDF6 library interface.
          Concrete implementation of the PDFModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  March 20, 2018

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENIE_LHAPDF6_INTERFACE_H_
#define _GENIE_LHAPDF6_INTERFACE_H_

#include "Physics/PartonDistributions/PDFModelI.h"

namespace LHAPDF
{
  class PDF;
}

namespace genie {

class LHAPDF6 : public PDFModelI {

public:

  LHAPDF6();
  LHAPDF6(string config);
  virtual ~LHAPDF6();

  // Implement the PDFModelI interface

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

  // Override the default "Configure" implementation 
  // of the Algorithm interface

  void Configure (const Registry & config);
  void Configure (string config);

private:

  void LoadConfig (void);

  string fSetName;
  int    fMemberID;

  LHAPDF::PDF * fLHAPDF;

};

}         // genie namespace
#endif    // _GENIE_LHAPDF6_INTERFACE_H_
