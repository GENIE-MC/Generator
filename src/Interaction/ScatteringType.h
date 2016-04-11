//____________________________________________________________________________
/*!

\class    genie::ScatteringType

\brief    Enumeration of scattering types

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 06, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________

#ifndef _SCATTERING_TYPE_H_
#define _SCATTERING_TYPE_H_

#include <cassert>
#include <string>

using std::string;

namespace genie {

typedef enum EScatteringType {

  kScNull = 0,
  kScQuasiElastic,
  kScSingleKaon,
  kScDeepInelastic,
  kScResonant,
  kScCoherent,
  kScDiffractive,
  kScNuElectronElastic,
  kScInverseMuDecay,
  kScAMNuGamma,
  kScMEC,
  kScCoherentElas,
  kScInverseBetaDecay,
  kScGlashowResonance,
  kScIMDAnnihilation

} ScatteringType_t;

class ScatteringType
{
public:

  //__________________________________________________________________________
  static string AsString(ScatteringType_t type)
  {
    switch (type) {

      case(kScQuasiElastic) :      return "QES";       break;
      case(kScDeepInelastic) :     return "DIS";       break;
      case(kScResonant) :          return "RES";       break;
      case(kScSingleKaon) :        return "1Kaon";     break;
      case(kScCoherent) :          return "COH";       break;
      case(kScDiffractive) :       return "DFR";       break;
      case(kScNuElectronElastic) : return "NuEEL";     break;
      case(kScInverseMuDecay) :    return "IMD";       break;
      case(kScAMNuGamma) :         return "AMNuGamma"; break;
      case(kScMEC) :               return "MEC";       break;
      case(kScCoherentElas) :      return "COHEl";     break;
      case(kScInverseBetaDecay) :  return "IBD";       break;
      case(kScGlashowResonance) :  return "GLR";       break;
      case(kScIMDAnnihilation) :   return "IMDAnh";    break;
      default :                    return "Unknown";   break;
    }
    return "Unknown";
  }
  //__________________________________________________________________________
};

}      // genie namespace
#endif // _SCATTERING_TYPE_H_
