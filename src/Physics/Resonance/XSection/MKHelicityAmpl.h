//____________________________________________________________________________
/*!

\class    genie::MKHelicityAmpl

\brief    A class holding the Rein-Sehgal's helicity amplitudes,
          extended for MK-model.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          MKHelicityAmplModelI interface, that it finds attached to itself.

\authors  Igor Kakorin <kakorin@jinr.ru> Joint Institute for Nuclear Research \n
          based on code by 
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Nov 12, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MK_HELICITY_AMPL_H_
#define _MK_HELICITY_AMPL_H_

#include <iostream>

#include <TMath.h>

#include "Framework/Interaction/Interaction.h"

using std::ostream;

namespace genie {

class MKHelicityAmpl;
ostream & operator<< (ostream & stream, const MKHelicityAmpl & hamp);

class MKHelicityAmpl {

friend class MKHelicityAmplModelCC;
friend class MKHelicityAmplModelNCp;
friend class MKHelicityAmplModelNCn;


public:

  MKHelicityAmpl();
  MKHelicityAmpl(const MKHelicityAmpl & hamp);
  ~MKHelicityAmpl() { }
  


  //! return helicity amplitude
  double AmpVMinus1  (void) const  { return fVMinus1;  } /* fV(-1)  */
  double AmpVPlus1   (void) const  { return fVPlus1;   } /* fV(+1)  */
  double AmpVMinus3  (void) const  { return fVMinus3;  } /* fV(-3)  */
  double AmpVPlus3   (void) const  { return fVPlus3;   } /* fV(+3)  */
  double AmpV0LMinus (void) const  { return fV0LMinus; } /* fV(0L-) */
  double AmpV0LPlus  (void) const  { return fV0LPlus;  } /* fV(0L+) */
  double AmpV0RMinus (void) const  { return fV0RMinus; } /* fV(0R-) */
  double AmpV0RPlus  (void) const  { return fV0RPlus;  } /* fV(0R+) */
  
  double AmpAMinus1  (void) const  { return fAMinus1;  } /* fA(-1)  */
  double AmpAPlus1   (void) const  { return fAPlus1;   } /* fA(+1)  */
  double AmpAMinus3  (void) const  { return fAMinus3;  } /* fA(-3)  */
  double AmpAPlus3   (void) const  { return fAPlus3;   } /* fA(+3)  */
  double AmpA0LMinus (void) const  { return fA0LMinus; } /* fA(0L-) */
  double AmpA0LPlus  (void) const  { return fA0LPlus;  } /* fA(0L+) */
  double AmpA0RMinus (void) const  { return fA0RMinus; } /* fA(0R-) */
  double AmpA0RPlus  (void) const  { return fA0RPlus;  } /* fA(0R+) */
  

  friend ostream & operator<< (ostream & stream, const MKHelicityAmpl & hamp);

  void Print(ostream & stream) const;

private:

  void   Init(void);

  double fVMinus1;
  double fVPlus1; 
  double fVMinus3;
  double fVPlus3; 
  double fV0LMinus;
  double fV0LPlus;
  double fV0RMinus;
  double fV0RPlus;

  double fAMinus1;
  double fAPlus1; 
  double fAMinus3;
  double fAPlus3; 
  double fA0LMinus;
  double fA0LPlus;
  double fA0RMinus;
  double fA0RPlus;
  
  
};

}        // genie namespace

#endif   // _MK_HELICITY_AMPL_H_


