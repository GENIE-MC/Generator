//____________________________________________________________________________
/*!

\class   genie::GVldContext

\brief   Validity Context for an Event Generator

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 20, 2004

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _GENERATOR_VALIDITY_CONTEXT_H_
#define _GENERATOR_VALIDITY_CONTEXT_H_

#include <string>
#include <iostream>

#include "Interaction/ScatteringType.h"
#include "Interaction/InteractionType.h"

using std::string;
using std::ostream;

namespace genie {

class Interaction;

class GVldContext {

public :
  GVldContext();
  GVldContext(const GVldContext & validity_context);
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
