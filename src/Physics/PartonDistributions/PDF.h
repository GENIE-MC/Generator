//____________________________________________________________________________
/*!

\class    genie::PDF

\brief    A class to store PDFs.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given (x,q^2) pair,
          that it then delegates to the algorithmic object, implementing the
          PDFModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PDF_H_
#define _PDF_H_

#include <iostream>

#include "Physics/PartonDistributions/PDFModelI.h"

using std::ostream;

namespace genie {

class PDF;
ostream & operator << (ostream & stream, const PDF & pdf_set);

class PDF {

public:

  PDF();
  PDF(const PDF & pdf_set);
  virtual ~PDF();

  //-- methods to set a PDFModelI and compute PDFs
  void   SetModel  (const PDFModelI * model);
  void   Calculate (double x, double q2);

  //-- methods to access the computed PDFs
  double UpValence   (void) const { return fUpValence;   }
  double DownValence (void) const { return fDownValence; }
  double UpSea       (void) const { return fUpSea;       }
  double DownSea     (void) const { return fDownSea;     }
  double Strange     (void) const { return fStrange;     }
  double Charm       (void) const { return fCharm;       }
  double Bottom      (void) const { return fBottom;      }
  double Top         (void) const { return fTop;         }
  double Gluon       (void) const { return fGluon;       }

  //-- methods to scale sea and valence PDFs (eg used to apply 
  //   corrections from non-QCD based fits / etc see Bodek Yang model)
  void ScaleValence     (double kscale);
  void ScaleSea         (double kscale);
  void ScaleUpValence   (double kscale);
  void ScaleDownValence (double kscale);
  void ScaleUpSea       (double kscale);
  void ScaleDownSea     (double kscale);
  void ScaleStrange     (double kscale);
  void ScaleCharm       (double kscale);

  //-- reseting/copying methods
  void Reset (void);
  void Copy  (const PDF & pdf_set);

  //-- printing methods & operators
  void Print(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const PDF & pdf_set);

protected:

  void Init(void);

  double fUpValence;
  double fDownValence;
  double fUpSea;
  double fDownSea;
  double fStrange;
  double fCharm;
  double fBottom;
  double fTop;
  double fGluon;

  const PDFModelI * fModel;
};

}         // genie namespace

#endif    // _PDF_H_
