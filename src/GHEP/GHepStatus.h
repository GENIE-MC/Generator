//____________________________________________________________________________
/*!

\class    genie::GHepStatus

\brief    GHepParticle Status

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 20, 2004

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
   kIstIntermediateState          =  2,
   kIstDecayedState               =  3,
   kIstNucleonTarget              = 11,
   kIstDISPreFragmHadronicState   = 12,
   kIstPreDecayResonantState      = 13,
   kIstHadronInTheNucleus         = 14

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
     case kIstIntermediateState:       
                            return  "<Intermediate State>"; break;
     case kIstDecayedState:            
                            return  "<Decayed State>"; break;
     case kIstNucleonTarget:            
                            return  "<Nucleon Target>"; break;
     case kIstDISPreFragmHadronicState: 
                            return  "<DIS Pre-Fragm. Hadronic State>"; break;
     case kIstPreDecayResonantState:   
                            return  "<Resonant Pre-Decayed State>"; break;
     case kIstHadronInTheNucleus:     
                            return  "<Hadron in the Nucleus>"; break;

     default:  break;
     }

     return "<->";
  }

};

}         // genie

#endif    // _STDHEP_STATUS_H_
