//____________________________________________________________________________
/*!

\class    genie::FKR

\brief    Simple struct-like class holding the Feynmann-Kislinger-Ravndall 
          (FKR) baryon excitation model parameters.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _FKR_H_
#define _FKR_H_

#include <iostream>

using std::ostream;

namespace genie {

class FKR {

public:

  friend ostream & operator<< (ostream & stream, const FKR & parameters);

  double Lamda;
  double Tv;
  double Rv;
  double S;
  double Ta;
  double Ra;
  double B;
  double C;
  double R;
  double T;
  double Tplus;
  double Tminus;
  double Rplus;
  double Rminus;

  void Reset (void);
  void Print (ostream & stream) const;

  FKR();
  ~FKR();
};

}        // genie namespace

#endif   // _FKR_H_
