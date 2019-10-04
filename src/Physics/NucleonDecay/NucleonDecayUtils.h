//____________________________________________________________________________
/*!

\class    genie::utils::nucleon_decay

\brief    Utilities for simulating nucleon decay

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 03, 2011

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NUCLEON_DECAY_UTILS_H_
#define _NUCLEON_DECAY_UTILS_H_

#include <string>

#include "Framework/GHEP/GHepStatus.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Physics/NucleonDecay/NucleonDecayMode.h"

using std::string;

namespace genie {
  namespace utils {
    namespace nucleon_decay {
      
      string         AsString              (NucleonDecayMode_t ndm, int npdg = 0);
      bool           IsValidMode           (NucleonDecayMode_t ndm, int npdg = 0);

      // The DecayedNucleonPdgCode utiliity is used ONLY for decay modes that do NOT require specifying a decayed nucleon PDG on the command-line. For these modes, this method returns that PDG code. Otherwise, the decayed nucleon PDG given on the command-line is used. So, no 2nd argument to this
      int            DecayedNucleonPdgCode (NucleonDecayMode_t ndm);

      PDGCodeList    DecayProductList      (NucleonDecayMode_t ndm, int npdg = 0);
      GHepStatus_t   DecayProductStatus    (bool in_nucleus, int pdgc);
      
    } // nucleon_decay namespace
  } // utils namespace
} // genie namespace

#endif // _NUCLEON_DECAY_UTILS_H_
