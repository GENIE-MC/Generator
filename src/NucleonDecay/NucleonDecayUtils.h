//____________________________________________________________________________
/*!

\class    genie::utils::nucleon_decay

\brief    Utilities for simulating nucleon decay

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 03, 2011

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NUCLEON_DECAY_UTILS_H_
#define _NUCLEON_DECAY_UTILS_H_

#include <string>

#include "GHEP/GHepStatus.h"
#include "NucleonDecay/NucleonDecayMode.h"
#include "PDG/PDGCodeList.h"

using std::string;

namespace genie {
 namespace utils {
  namespace nucleon_decay {

     string         AsString              (NucleonDecayMode_t ndm);
     bool           IsValidMode           (NucleonDecayMode_t ndm);
     int            DecayedNucleonPdgCode (NucleonDecayMode_t ndm);
     PDGCodeList    DecayProductList      (NucleonDecayMode_t ndm);
     GHepStatus_t   DecayProductStatus    (bool in_nucleus, int pdgc);

  } // nucleon_decay namespace
 } // utils namespace
} // genie namespace

#endif // _NUCLEON_DECAY_UTILS_H_
