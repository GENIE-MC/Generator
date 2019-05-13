//____________________________________________________________________________
/*!

\class    genie::utils::neutron_osc

\brief    Utilities for simulating neutron oscillation

\author   Jeremy Hewes, Georgia Karagiorgi
          University of Manchester

          Adapted from the NucleonDecay package (Author: Costas Andreopoulos).

\created  November, 2016

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NNBAR_OSC_UTILS_H_
#define _NNBAR_OSC_UTILS_H_

#include <string>

#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Physics/NNBarOscillation/NNBarOscMode.h"

using std::string;

namespace genie {
 namespace utils {
  namespace nnbar_osc {

     string         AsString                   (NNBarOscMode_t ndm);
     bool           IsValidMode                (NNBarOscMode_t ndm);
     int            AnnihilatingNucleonPdgCode (NNBarOscMode_t ndm);
     PDGCodeList    DecayProductList           (NNBarOscMode_t ndm);
     GHepStatus_t   DecayProductStatus         (bool in_nucleus, int pdgc);

  } // nnbar_osc namespace
 } // utils namespace
} // genie namespace

#endif // _NNBAR_OSC_UTILS_H_
