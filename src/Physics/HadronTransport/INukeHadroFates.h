//____________________________________________________________________________
/*!

\class    genie::INukeHadroFates

\brief    An enumeration of possible hadron "fates" taken into account by the
          INTRANUKE hadron transport MC.

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
          Costas Andreopoulos <c.andreopoulos \at cern.ch>, Rutherford Lab.

\created  November 1, 2005

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
           
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_FATES_H_
#define _INTRANUKE_FATES_H_

#include <string>

namespace genie {

using std::string;

// Fates in legacy INTRANUKE's HN mode
//
typedef enum EINukeFateHNLeg_t {

   kIHNLegFtUndefined = 0, 
   kIHNLegFtNoInteraction,
   kIHNLegFtCEx,       // cex
   kIHNLegFtElas,      // elas
   kIHNLegFtInelas,    // inelas
   kIHNLegFtAbs,       // abs 
   kIHNLegFtCmp         //cmp

} INukeFateHNLeg_t;   

// Fates in legacy INTRANUKE's HA mode
//
typedef enum EINukeFateHALeg_t {

   kIHALegFtUndefined = 0,
   kIHALegFtNoInteraction,  // no interaction 
   kIHALegFtCEx,            // cex
   kIHALegFtElas,           // elas
   kIHALegFtInelas,         // inelas
   kIHALegFtAbs,            // abs
   kIHALegFtKo, 	         // knock out
   kIHALegFtCmp,            // compound nucleus
   kIHALegFtPiProd,         // pi production
   kIHALegFtInclPip,        // pi production : inclusive pi+
   kIHALegFtInclPim,        // pi production : inclusive pi-
   kIHALegFtInclPi0,        // pi production : inclusive pi0 
   kIHALegFtDCEx            // dcex

} INukeFateHALeg_t;   

class INukeHadroFates {

public:
  //__________________________________________________________________________
  static string AsString(INukeFateHNLeg_t fate) {
     switch (fate) {
      case kIHNLegFtUndefined : return "** Undefined HN-mode fate **"; break;
      case kIHNLegFtCEx       : return "HN-mode / cex";    break;
      case kIHNLegFtElas      : return "HN-mode / elas";   break;
      case kIHNLegFtInelas    : return "HN-mode / inelas"; break;
      case kIHNLegFtAbs       : return "HN-mode / abs";    break;
      case kIHNLegFtCmp	   : return "HN-mode / compound"; break;
      case kIHNLegFtNoInteraction : return "HN-mode / no interaction"; break;
      default              : break; 
     }
     return "** Undefined HN-mode fate **"; 
  }
  //__________________________________________________________________________
  static string AsString(INukeFateHALeg_t fate) {
     switch (fate) {
      case kIHALegFtUndefined : return "** Undefined HA-mode fate **"; break;
      case kIHALegFtNoInteraction : return "HA-mode / no interaction"; break;
      case kIHALegFtCEx       : return "HA-mode / cex";            break;
      case kIHALegFtElas      : return "HA-mode / elas";           break;
      case kIHALegFtInelas    : return "HA-mode / inelas";         break;
      case kIHALegFtAbs       : return "HA-mode / abs";            break;
      case kIHALegFtKo        : return "HA-mode / knock-out";      break;
      case kIHALegFtCmp       : return "HA-mode / compound";       break;
      case kIHALegFtPiProd    : return "HA-mode / pi-production" ; break;
      case kIHALegFtInclPip   : return "HA-mode / pi-prod incl pi+";   break;
      case kIHALegFtInclPim   : return "HA-mode / pi-prod incl pi-";   break;
      case kIHALegFtInclPi0   : return "HA-mode / pi-prod incl pi0";   break;
      case kIHALegFtDCEx      : return "HA-mode / dcex";           break;
      default              : break;
     }
     return "** Undefined HA-mode fate **"; 
  }
  //__________________________________________________________________________
  static string AsSimpleString(INukeFateHALeg_t fate) {
     switch (fate) {
      case kIHALegFtUndefined : return "undefined"; break;
      case kIHALegFtNoInteraction : return "no interaction"; break;
      case kIHALegFtCEx       : return "cex";            break;
      case kIHALegFtElas      : return "elas";           break;
      case kIHALegFtInelas    : return "inelas";         break;
      case kIHALegFtAbs       : return "abs";            break;
      case kIHALegFtKo        : return "knock out";      break; 
      case kIHALegFtCmp       : return "compound";       break;
      case kIHALegFtPiProd    : return "pi prod";          break;
      case kIHALegFtDCEx      : return "dcex";             break;
      default              : break;
     }
     return "undefined"; 
  }
  //__________________________________________________________________________

};

}      // genie
#endif // _INTRANUKE_FATES_H_
