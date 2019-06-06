//____________________________________________________________________________
/*!

\class    genie::INukeHadroFates

\brief    An enumeration of possible hadron "fates" taken into account by the
          INTRANUKE hadron transport MC.

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>, Rutherford Lab.

\created  November 1, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_FATES_H_
#define _INTRANUKE_FATES_H_

#include <string>



namespace genie {

using std::string;

// Fates in INTRANUKE's HN mode
//
typedef enum EINukeFateHN_t {

   kIHNFtUndefined = 0, 
   kIHNFtNoInteraction,
   kIHNFtCEx,       // cex
   kIHNFtElas,      // elas
   kIHNFtInelas,    // inelas
   kIHNFtAbs,       // abs 
   kIHNFtCmp         //cmp

} INukeFateHN_t;   

// Fates in INTRANUKE's HA mode
//
typedef enum EINukeFateHA_t {

   kIHAFtUndefined = 0,
   kIHAFtNoInteraction,  // no interaction 
   kIHAFtCEx,            // cex
   kIHAFtElas,           // elas
   kIHAFtInelas,         // inelas
   kIHAFtAbs,            // abs
   kIHAFtKo, 	         // knock out
   kIHAFtCmp,            // compound nucleus
   kIHAFtPiProd,         // pi production
   kIHAFtInclPip,        // pi production : inclusive pi+
   kIHAFtInclPim,        // pi production : inclusive pi-
   kIHAFtInclPi0,        // pi production : inclusive pi0 
   kIHAFtDCEx            // dcex

} INukeFateHA_t;   

class INukeHadroFates {

public:
  //__________________________________________________________________________
  static string AsString(INukeFateHN_t fate) {
     switch (fate) {
      case kIHNFtUndefined : return "** Undefined HN-mode fate **"; break;
      case kIHNFtCEx       : return "HN-mode / cex";    break;
      case kIHNFtElas      : return "HN-mode / elas";   break;
      case kIHNFtInelas    : return "HN-mode / inelas"; break;
      case kIHNFtAbs       : return "HN-mode / abs";    break;
      case kIHNFtCmp	   : return "HN-mode / compound"; break;
      case kIHNFtNoInteraction : return "HN-mode / no interaction"; break;
      default              : break; 
     }
     return "** Undefined HN-mode fate **"; 
  }
  //__________________________________________________________________________
  static string AsString(INukeFateHA_t fate) {
     switch (fate) {
      case kIHAFtUndefined : return "** Undefined HA-mode fate **"; break;
      case kIHAFtNoInteraction : return "HA-mode / no interaction"; break;
      case kIHAFtCEx       : return "HA-mode / cex";            break;
      case kIHAFtElas      : return "HA-mode / elas";           break;
      case kIHAFtInelas    : return "HA-mode / inelas";         break;
      case kIHAFtAbs       : return "HA-mode / abs";            break;
      case kIHAFtKo        : return "HA-mode / knock-out";      break;
      case kIHAFtCmp       : return "HA-mode / compound";       break;
      case kIHAFtPiProd    : return "HA-mode / pi-production" ; break;
      case kIHAFtInclPip   : return "HA-mode / pi-prod incl pi+";   break;
      case kIHAFtInclPim   : return "HA-mode / pi-prod incl pi-";   break;
      case kIHAFtInclPi0   : return "HA-mode / pi-prod incl pi0";   break;
      case kIHAFtDCEx      : return "HA-mode / dcex";           break;
      default              : break;
     }
     return "** Undefined HA-mode fate **"; 
  }
  //__________________________________________________________________________
  static string AsSimpleString(INukeFateHA_t fate) {
     switch (fate) {
      case kIHAFtUndefined : return "undefined"; break;
      case kIHAFtNoInteraction : return "no interaction"; break;
      case kIHAFtCEx       : return "cex";            break;
      case kIHAFtElas      : return "elas";           break;
      case kIHAFtInelas    : return "inelas";         break;
      case kIHAFtAbs       : return "abs";            break;
      case kIHAFtKo        : return "knock out";      break; 
      case kIHAFtCmp       : return "compound";       break;
      case kIHAFtPiProd    : return "pi prod";          break;
      case kIHAFtDCEx      : return "dcex";             break;
      default              : break;
     }
     return "undefined"; 
  }
  //__________________________________________________________________________
};

}      // genie
#endif // _INTRANUKE_FATES_H_
