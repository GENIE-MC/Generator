//____________________________________________________________________________
/*!

\class    genie::ScatteringType

\brief    Enumeration of scattering types : QEL, DIS, RES, ...

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004
 
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
  kScElastic,
  kScQuasiElastic,
  kScDeepInelastic,
  kScResonant,
  kScCoherent,
  kScInverseMuDecay

} ScatteringType_t;

class ScatteringType
{
public:

  //__________________________________________________________________________
  static string AsString(ScatteringType_t type)
  {
    switch (type) {

      case(kScElastic) :        return "ELS";      break;
      case(kScQuasiElastic) :   return "QES";      break;
      case(kScDeepInelastic) :  return "DIS";      break;
      case(kScResonant) :       return "RES";      break;
      case(kScCoherent) :       return "COH";      break;
      case(kScInverseMuDecay) : return "IMD";      break;
      default :                 return "Unknown";  break;
    }
    return "Unknown";
  }
  //__________________________________________________________________________
  static ScatteringType_t FromString(string type)
  {
    //-- Make uppercase/lowercase irrelevant

    for(unsigned int i=0; i<type.size(); i++) type[i] = toupper(type[i]);

    //-- Figure out the ScatteringType_t from the input string

    const char * t = type.c_str();

    if ( strcmp(t,"ELASTIC") == 0 || strcmp(t,"ELS") == 0 ) return kScElastic;

    else if ( strcmp(t,"QUASIELASTIC") == 0 ||
                    strcmp(t,"QUASI-ELASTIC") == 0 ||
                          strcmp(t,"QUASI ELASTIC") == 0 ||
                                strcmp(t,"QEL") == 0 ) return kScQuasiElastic;

    else if ( strcmp(t,"DEEPINELASTIC") == 0 ||
                  strcmp(t,"DEEP-INELASTIC") == 0 ||
                       strcmp(t,"DEEP INELASTIC") == 0 ||
                               strcmp(t,"DIS") == 0 ) return kScDeepInelastic;

    else if ( strcmp(t,"RESONANCE") == 0 ||
                  strcmp(t,"SINGLE PION PRODUCTION") == 0 ||
                       strcmp(t,"SINGLE-PION-PRODUCTION") == 0 ||
                           strcmp(t,"SPP") == 0 ||
                                    strcmp(t,"RES") == 0 ) return kScResonant;

    else if ( strcmp(t,"COHERENT") == 0 ||
                   strcmp(t,"COHERENT PION PRODUCTION") == 0 ||
                         strcmp(t,"COHERENT-PION-PRODUCTION") == 0 ||
                                    strcmp(t,"COH") == 0 ) return kScCoherent;

    else if ( strcmp(t,"INVERSE-MUON-DECAY") == 0 ||
                   strcmp(t,"INVERSE MUON DECAY") == 0 ||
                              strcmp(t,"IMD") == 0 ) return kScInverseMuDecay;

    else return kScNull;
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _SCATTERING_TYPE_H_
