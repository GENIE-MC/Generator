//____________________________________________________________________________
/*!

\class    genie::utils::neutron_osc

\brief    Utilities for simulating neutron oscillation

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 03, 2011

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NEUTRON_OSC_UTILS_H_
#define _NEUTRON_OSC_UTILS_H_

#include <string>

#include "GHEP/GHepStatus.h"
#include "NeutronOsc/NeutronOscMode.h"
#include "PDG/PDGCodeList.h"

using std::string;

namespace genie {
 namespace utils {
  namespace neutron_osc {

     string         AsString                   (NeutronOscMode_t ndm);
     bool           IsValidMode                (NeutronOscMode_t ndm);
     int            AnnihilatingNucleonPdgCode (NeutronOscMode_t ndm);
     PDGCodeList    DecayProductList           (NeutronOscMode_t ndm);
     GHepStatus_t   DecayProductStatus         (bool in_nucleus, int pdgc);

  } // neutron_osc namespace
 } // utils namespace
} // genie namespace

#endif // _NEUTRON_OSC_UTILS_H_
