//____________________________________________________________________________
/*!

\class    genie::RESVectFFAmplitude

\brief    A class holding the RES Vector FF Amplitudes

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  January 2023

\cpright  Copyright (c) 2023-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _RES_VECT_FF_AMPL_H_
#define _RES_VECT_FF_AMPL_H_

#include <iostream>

#include <TMath.h>

using std::ostream;

namespace genie {

class RESVectFFAmplitude {

public:

  RESVectFFAmplitude();
  RESVectFFAmplitude(const RESVectFFAmplitude & hamp);
  ~RESVectFFAmplitude() { }

  void SetAmplA12 ( const double a12 ) { fA12 = a12 ; }
  void SetAmplA32 ( const double a32 ) { fA32 = a32; }
  void SetAmplS12 ( const double s12 ) { fS12 = s12 ; }

  //! return helicity amplitude
  double AmplA12 (void) const { return fA12 ; }
  double AmplA32 (void) const { return fA32 ; }
  double AmplS12 (void) const { return fS12 ; }

  double Ampl2A12 (void) const { return pow(fA12,2) ; }
  double Ampl2A32 (void) const { return pow(fA32,2) ; }
  double Ampl2S12 (void) const { return pow(fS12,2) ; }

  friend ostream & operator<< (ostream & stream, const RESVectFFAmplitude & hamp);

  void Print(ostream & stream) const;

private:

  void   Init(void);

  double fA12 ; 
  double fA32 ; 
  double fS12 ; 
};

ostream & operator<< (ostream & stream, const RESVectFFAmplitude & hamp);

}        // genie namespace

#endif   // _RES_VECT_FF_AMPL_H_
