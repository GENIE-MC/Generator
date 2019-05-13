//____________________________________________________________________________
/*!

\class    genie::ScatteringType

\brief    Enumeration of scattering types

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Changes required to implement the GENIE Boosted Dark Matter module
          were installed by Josh Berger (Univ. of Wisconsin)

\created  May 06, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
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

// Note: please attach new _neutrino_ scattering modes to the _end_ of the
// list of neutrino enums, and new dark matter modes to the end of the list
// of dark matter enums, etc. If adding an entirely new set of enums, please
// append to the end of the total list and set a new enum counter value.
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
  kScIMDAnnihilation,
  kScDarkMatterElastic = 101,
  kScDarkMatterDeepInelastic

} ScatteringType_t;

class ScatteringType
{
public:

  //__________________________________________________________________________
  static string AsString(ScatteringType_t type)
  {
    switch (type) {

      case(kScQuasiElastic) :            return "QES";       break;
      case(kScSingleKaon) :              return "1Kaon";     break;
      case(kScDeepInelastic) :           return "DIS";       break;
      case(kScResonant) :                return "RES";       break;
      case(kScCoherent) :                return "COH";       break;
      case(kScDiffractive) :             return "DFR";       break;
      case(kScNuElectronElastic) :       return "NuEEL";     break;
      case(kScInverseMuDecay) :          return "IMD";       break;
      case(kScAMNuGamma) :               return "AMNuGamma"; break;
      case(kScMEC) :                     return "MEC";       break;
      case(kScCoherentElas) :            return "COHEl";     break;
      case(kScInverseBetaDecay) :        return "IBD";       break;
      case(kScGlashowResonance) :        return "GLR";       break;
      case(kScIMDAnnihilation) :         return "IMDAnh";    break;
      case(kScDarkMatterElastic) :       return "DME";       break;
      case(kScDarkMatterDeepInelastic) : return "DMDIS";     break;
      default :                          return "Unknown";   break;
    }
    return "Unknown";
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _SCATTERING_TYPE_H_
