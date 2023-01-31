//____________________________________________________________________________
/*!

\class    genie::MAIDHelicityAmpl

\brief    A class holding the Rein-Sehgal's helicity amplitudes.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          MAIDHelicityAmplModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _MAID_HELICITY_AMPL_H_
#define _MAID_HELICITY_AMPL_H_

#include <iostream>

#include <TMath.h>

using std::ostream;

namespace genie {

class MAIDHelicityAmpl;
ostream & operator<< (ostream & stream, const MAIDHelicityAmpl & hamp);

class MAIDHelicityAmpl {

friend class MAIDHelicityAmplModelEMp;
friend class MAIDHelicityAmplModelEMn;

public:

  MAIDHelicityAmpl();
  MAIDHelicityAmpl(const MAIDHelicityAmpl & hamp);
  ~MAIDHelicityAmpl() { }

  //! return helicity amplitude
  double AmplA12 (void) const { return fA12 ; }
  double AmplA32 (void) const { return fA32 ; }
  double AmplS12 (void) const { return fS12 ; }

  friend ostream & operator<< (ostream & stream, const MAIDHelicityAmpl & hamp);

  void Print(ostream & stream) const;

private:

  void   Init(void);

  double fA12 ; 
  double fA32 ; 
  double fS12 ; 
};

}        // genie namespace

#endif   // _MAID_HELICITY_AMPL_H_
