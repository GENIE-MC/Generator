//____________________________________________________________________________
/*!

\class    genie::FKR

\brief    Simple struct-like class holding the Feynmann-Kislinger-Ravndall
          (FKR) baryon excitation model parameters.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 03, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _FKR_H_
#define _FKR_H_

#include <iostream>

using std::ostream;

namespace genie {

class FKR;
ostream & operator<< (ostream & stream, const FKR & parameters);

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
