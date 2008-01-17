//____________________________________________________________________________
/*!

\class    genie::ScatteringType

\brief    Enumeration of scattering types

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 06, 2004

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
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
  kScDeepInelastic,
  kScResonant,
  kScCoherentPiProd,
  kScCoherentElas,
  kScDiffractive,
  kScNuElectronElastic,
  kScInverseMuDecay

} ScatteringType_t;

class ScatteringType
{
public:

  //__________________________________________________________________________
  static string AsString(ScatteringType_t type)
  {
    switch (type) {

      case(kScQuasiElastic) :      return "QES";      break;
      case(kScDeepInelastic) :     return "DIS";      break;
      case(kScResonant) :          return "RES";      break;
      case(kScCoherentPiProd) :    return "COHPi";    break;
      case(kScCoherentElas) :      return "COHEl";    break;
      case(kScDiffractive) :       return "DFR";      break;
      case(kScNuElectronElastic) : return "NuEEL";    break;
      case(kScInverseMuDecay) :    return "IMD";      break;
      default :                    return "Unknown";  break;
    }
    return "Unknown";
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _SCATTERING_TYPE_H_
