//____________________________________________________________________________
/*!

\class    genie::PDF

\brief    A class to store PDFs.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given (x,q^2) pair,
          that it then delegates to the algorithmic object, implementing the
          PDFModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#ifndef _PDF_H_
#define _PDF_H_

#include <iostream>

#include "PDF/PDFModelI.h"

using std::ostream;

namespace genie {

class PDF {

public:

  PDF();
  PDF(const PDF & pdf_set);
  virtual ~PDF();

  virtual void   SetModel  (const PDFModelI * model);
  virtual void   Calculate (double x, double q2);

  double UpValence   (void) const { return fUpValence;   }
  double DownValence (void) const { return fDownValence; }
  double UpSea       (void) const { return fUpSea;       }
  double DownSea     (void) const { return fDownSea;     }
  double Strange     (void) const { return fStrange;     }
  double Charm       (void) const { return fCharm;       }
  double Bottom      (void) const { return fBottom;      }
  double Top         (void) const { return fTop;         }
  double Gluon       (void) const { return fGluon;       }

  void Print(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const PDF & pdf_set);

protected:

  void InitPDFs(void);

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
