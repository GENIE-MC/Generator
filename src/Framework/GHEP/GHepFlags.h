//____________________________________________________________________________
/*!

\class    genie::GHepFlag

\brief    An enumeration of event flags. Each represents a physical condition 
          or a computational error. If any is set the event would be marked as 
          unphysical.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 06, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GHEP_FLAGS_H_
#define _GHEP_FLAGS_H_

namespace genie {

  typedef enum EGHepFlag {

     kGenericErr      = 0,
     kPauliBlock      = 1,
     kBelowThrNRF     = 2,
     kBelowThrERF     = 3,
     kKineGenErr      = 4,
     kHadroSysGenErr  = 5,
     kLeptoGenErr     = 6,
     kDecayErr        = 7

  } GHepFlag_t;

class GHepFlags {

 public:
  //__________________________________________________________________________
  static const char * Describe(GHepFlag_t flag) 
  {
     switch (flag) {
     case kGenericErr :
            return "Generic error";
            break;
     case kPauliBlock :
            return "Pauli-blocked event";
            break;
     case kBelowThrNRF :
            return "E<Ethr in hit nucleon rest frame";
            break;
     case kBelowThrERF :
            return "E<Ethr in hit e- rest frame";
            break;
     case kKineGenErr :
            return "Generic error in kinematic generation";
            break;
     case kHadroSysGenErr : 
            return "Generic error in f/s hadronic system generation";
            break;
     case kLeptoGenErr : 
            return "Generic error in f/s lepton generation";
            break;
     case kDecayErr : 
            return "Generic error during unstable particle decay";
            break;
     default:
            return "Unknown GHEP flag";
            break;
     }
     return "Unknown GHEP flag";
  }
  //__________________________________________________________________________
  static unsigned int NFlags(void) { return 16; }
  //__________________________________________________________________________
};

}        // genie namespace

#endif   // _GHEP_FLAGS_H_
