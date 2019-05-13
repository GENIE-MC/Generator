//____________________________________________________________________________
/*!

\class    genie::mueloss::MuELProcess

\brief    Enumeration of muon energy loss processes

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  December 10, 2003

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MUELOSS_PROCESS_H_
#define _MUELOSS_PROCESS_H_

#include "Framework/Conventions/Constants.h"

using namespace genie::constants;

namespace genie   {
namespace mueloss {

typedef enum EMuELProcess {

  eMupUndefined = 0,
  eMupIonization,
  eMupPairProduction,
  eMupBremsstrahlung,     
  eMupNuclearInteraction,
  eMupSum,

} MuELProcess_t;

class MuELProcess
{
public:
  //__________________________________________________________________________
  static const char * AsString(MuELProcess_t p) {
     switch(p) {
     case eMupIonization        : return "ionization";                  break;
     case eMupPairProduction    : return "direct e+e- pair production"; break;
     case eMupBremsstrahlung    : return "bremsstrahlung";              break;
     case eMupNuclearInteraction: return "nuclear interaction";         break;
     case eMupSum               : return "all";                         break;
     case eMupUndefined: 
     default:
         return "undefined process"; break;
     }
     return "undefined process";
  }
  //__________________________________________________________________________
  static double Threshold(MuELProcess_t p) {
     switch(p) {
     case eMupIonization        : return kMuonMass;                 break;
     case eMupPairProduction    : return kMuonMass+2*kElectronMass; break;
     case eMupBremsstrahlung    : return kMuonMass;                 break;
     case eMupNuclearInteraction: return kMuonMass;                 break;
     case eMupSum               : return kMuonMass;                 break;
     case eMupUndefined: 
     default:
         return 999999999; break;
     }
     return 99999999;
  }
  //__________________________________________________________________________
};

}      // mueloss namespace
}      // genie   namespace

#endif // _MUELOSS_PROCESS_H_
