//____________________________________________________________________________
/*!

\class    genie::GHepStatus

\brief    GHepParticle Status

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 20, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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
   kIStInitialState               =  0,   /* generator-level initial state */
   kIStStableFinalState           =  1,   /* generator-level final state: particles to be tracked by detector-level MC */
   kIStIntermediateState          =  2,
   kIStDecayedState               =  3,
   kIStCorrelatedNucleon          = 10,
   kIStNucleonTarget              = 11,
   kIStDISPreFragmHadronicState   = 12,
   kIStPreDecayResonantState      = 13,
   kIStHadronInTheNucleus         = 14,   /* hadrons inside the nucleus: marked for hadron transport modules to act on */
   kIStFinalStateNuclearRemnant   = 15,   /* low energy nuclear fragments entering the record collectively as a 'hadronic blob' pseudo-particle */
   kIStNucleonClusterTarget       = 16    // for composite nucleons before phase space decay
} 
GHepStatus_t; 
  
class GHepStatus {
 public:

  static const char * AsString(GHepStatus_t Ist) {
     switch (Ist) {
     case kIStUndefined:                
           return  "[undefined status]"; 
           break;
     case kIStInitialState:             
           return  "[initial state]"; 
           break;
     case kIStStableFinalState:         
           return  "[stable final state]"; 
           break;
     case kIStIntermediateState:       
           return  "[intermediate state]"; 
           break;
     case kIStDecayedState:            
           return  "[decayed state]"; 
           break;
     case kIStCorrelatedNucleon:
           return  "[other energetic initial state nucleons]"; 
           break;
     case kIStNucleonTarget:            
           return  "[nucleon target]"; 
           break;
     case kIStDISPreFragmHadronicState: 
           return  "[DIS pre-fragm. hadronic state]"; 
           break;
     case kIStPreDecayResonantState:   
           return  "[resonant pre-decayed state]"; 
           break;
     case kIStHadronInTheNucleus:     
           return  "[hadron in the nucleus]"; 
           break;
     case kIStFinalStateNuclearRemnant:     
           return  "[nuclear remnant]"; 
           break;
     case kIStNucleonClusterTarget:     
           return  "[nucleon cluster target]"; 
           break;
     default:  break;
     }
     return "[-]";
  }
};

}         // genie
#endif    // _STDHEP_STATUS_H_
