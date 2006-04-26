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
#include <vector>
#include <iostream>

#include "Interaction/ScatteringType.h"
#include "Interaction/InteractionType.h"

using std::string;
using std::vector;
using std::ostream;

namespace genie {

class Interaction;

class GVldContext {

public :

  GVldContext();
  GVldContext(const GVldContext & validity_context);
  ~GVldContext();

  void   Decode  ( string encoded_values );

  bool   IsValid (const Interaction * proc) const;
  bool   Clashes (const GVldContext & vld_context) const;  
 
  double Emin    (void) const { return fEmin; }
  double Emax    (void) const { return fEmax; }
  
  void   Print   (ostream & stream) const;
 
  friend ostream & operator<< (ostream & stream, const GVldContext & vldc);

private:

  void Init(void);

  void DecodePROC   ( string encoded_values );
  void DecodeCURR   ( string encoded_values );
  void DecodePROBE  ( string encoded_values );
  void DecodeTARGET ( string encoded_values );
  void DecodeENERGY ( string encoded_values );
    
  ScatteringType_t            fProc;     // process: QEL, DIS, RES,...
  vector<InteractionType_t> * fCurr;     // current: CC, NC, E/M
  double                      fEmin;     // min probe energy in validity range
  double                      fEmax;     // max probe energy in validity range
  vector<int> *               fProbes;   // PDG codes for valid probes
  vector<int> *               fTargets;  // PDG codes for valid targets
};

}      // genie namespace

#endif // _GENERATOR_VALIDITY_CONTEXT_H_
