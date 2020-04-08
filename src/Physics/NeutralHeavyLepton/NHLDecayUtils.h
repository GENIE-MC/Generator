//____________________________________________________________________________
/*!

\class    genie::utils::nhl

\brief    Utilities for simulating the decay of Neutral Heavy Leptons 

\author   
          
\created  November 03, 2011

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _NHL_DECAY_UTILS_H_
#define _NHL_DECAY_UTILS_H_

#include <string>

#include "Framework/ParticleData/PDGCodeList.h"
#include "Physics/NeutralHeavyLepton/NHLDecayMode.h"

using std::string;

namespace genie 
{
 namespace utils 
 {
  namespace nhl 
  {
      string       AsString                (NHLDecayMode_t nhldm);
      bool         IsKinematicallyAllowed  (NHLDecayMode_t nhldm, double Mnhl);
      PDGCodeList  DecayProductList        (NHLDecayMode_t nhldm);

   } // nhl 
 } // utils 
} // genie 

#endif // _NHL_DECAY_UTILS_H_
