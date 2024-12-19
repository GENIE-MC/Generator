//____________________________________________________________________________
/*!

\class   genie::GVldContext

\brief   Validity Context for an Event Generator

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created November 20, 2004

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _GENERATOR_VALIDITY_CONTEXT_H_
#define _GENERATOR_VALIDITY_CONTEXT_H_

#include <string>
#include <iostream>

#include "Framework/Interaction/ScatteringType.h"
#include "Framework/Interaction/InteractionType.h"

using std::string;
using std::ostream;

namespace genie {

class GVldContext;
class Interaction;

ostream & operator<< (ostream & stream, const GVldContext & vldc);

class GVldContext {

public :
  GVldContext();
 ~GVldContext();

  void   Decode  (string encoded_values);

  double Emin    (void) const { return fEmin; }
  double Emax    (void) const { return fEmax; }

  void   Print   (ostream & stream) const;

  friend ostream & operator<< (ostream & stream, const GVldContext & vldc);

private:

  void Init(void);

  void DecodeENERGY (string encoded_values);

  double fEmin;  // min probe energy in validity range
  double fEmax;  // max probe energy in validity range
};

}      // genie namespace

#endif // _GENERATOR_VALIDITY_CONTEXT_H_
