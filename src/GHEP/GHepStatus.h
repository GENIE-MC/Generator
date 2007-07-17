//____________________________________________________________________________
/*!

\class    genie::GHepStatus

\brief    GHepParticle Status

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  November 20, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/ 
//____________________________________________________________________________

#ifndef _STDHEP_STATUS_H_
#define _STDHEP_STATUS_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

typedef enum EGHepStatus {

   kIStUndefined                  = -1, 
   kIStInitialState               =  0,
   kIStStableFinalState           =  1,
   kIStIntermediateState          =  2,
   kIStDecayedState               =  3,
   kIStNucleonTarget              = 11,
   kIStDISPreFragmHadronicState   = 12,
   kIStPreDecayResonantState      = 13,
   kIStHadronInTheNucleus         = 14

} GHepStatus_t; 
  

class GHepStatus {

 public:

  static char * AsString(GHepStatus_t Ist) {

     switch (Ist) {

     case kIStUndefined:                
                            return  "<Undefined Status>"; break;
     case kIStInitialState:             
                            return  "<Initial State>"; break;
     case kIStStableFinalState:         
                            return  "<Stable Final State>"; break;
     case kIStIntermediateState:       
                            return  "<Intermediate State>"; break;
     case kIStDecayedState:            
                            return  "<Decayed State>"; break;
     case kIStNucleonTarget:            
                            return  "<Nucleon Target>"; break;
     case kIStDISPreFragmHadronicState: 
                            return  "<DIS Pre-Fragm. Hadronic State>"; break;
     case kIStPreDecayResonantState:   
                            return  "<Resonant Pre-Decayed State>"; break;
     case kIStHadronInTheNucleus:     
                            return  "<Hadron in the Nucleus>"; break;

     default:  break;
     }

     return "<->";
  }

};

}         // genie

#endif    // _STDHEP_STATUS_H_
