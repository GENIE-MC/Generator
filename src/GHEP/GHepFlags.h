//____________________________________________________________________________
/*!

\class    genie::GHepFlag

\brief    An enumeration of event flags. Each represents a physical condition 
          or a computational error. If any is set the event would be marked as 
          unphysical.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

*/
//____________________________________________________________________________

#ifndef _GHEP_FLAGS_H_
#define _GHEP_FLAGS_H_

namespace genie {

  typedef enum EGHepFlag {

     kGenericErr               = 0,
     kPauliBlock               = 1,
     kBelowThrNucleonRestFrame = 2,
     kBelowThrElecRestFrame    = 3,
     kNoAvailablePhaseSpace    = 4,
     kNoValidKinematics        = 5

  } GHepFlag_t;

class GHepFlags {

 public:
  //___________________________________________________
  static char * Describe(GHepFlag_t flag) 
  {
     switch (flag) {
     case kGenericErr :
            return "Generic error";
            break;
     case kPauliBlock :
            return "Pauli-blocking";
            break;
     case kBelowThrNucleonRestFrame :
            return "E<Ethr in N rest frame";
            break;
     case kBelowThrElecRestFrame :
            return "E<Ethr in e- rest frame";
            break;
     case kNoAvailablePhaseSpace :
            return "No available phase space";
            break;
     case kNoValidKinematics : 
            return "No kinematics selection";
            break;
     default:
            return "Unknown GHEP flag";
            break;
     }
  }
  //___________________________________________________
  static unsigned int NFlags(void) 
  {
     return 16;
  }
  //___________________________________________________
};

}        // genie namespace

#endif   // _GHEP_FLAGS_H_
