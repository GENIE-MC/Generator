//____________________________________________________________________________
/*!

\class    genie::INukeHadroFates

\brief    An enumeration of possible hadron "fates" taken into account by the
          INTRANUKE hadron transport MC.

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
          Costas Andreopoulos <c.andreopoulos \at cern.ch>, Rutherford Lab.

\created  September 22, 2025

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

     @ Sept, 2025 Copied by Mohamed Ismail <msi10@pitt.edu>, no changes
 
           
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_FATES_2025_H_
#define _INTRANUKE_FATES_2025_H_

#include <string>

using std::string;

namespace genie {

// Fates in INTRANUKE's HN mode
//
typedef enum EINukeFateHN2025_t {

   kIHN25FtUndefined = 0, 
   kIHN25FtNoInteraction,
   kIHN25FtCEx,       // cex
   kIHN25FtElas,      // elas
   kIHN25FtInelas,    // inelas
   kIHN25FtAbs,       // abs 
   kIHN25FtCmp         //cmp

} INukeFateHN2025_t;   

// Fates in INTRANUKE's HA mode
//
typedef enum EINukeFateHA2025_t {

   kIHA25FtUndefined = 0,
   kIHA25FtNoInteraction,  // no interaction 
   kIHA25FtCEx,            // cex
   //   kIHA25FtElas,           // elas
   kIHA25FtInelas,         // inelas
   kIHA25FtAbs,            // abs
   kIHA25FtKo, 	         // knock out
   kIHA25FtCmp,            // compound nucleus
   kIHA25FtPiProd,         // pi production
   kIHA25FtInclPip,        // pi production : inclusive pi+
   kIHA25FtInclPim,        // pi production : inclusive pi-
   kIHA25FtInclPi0,        // pi production : inclusive pi0 
   kIHA25FtDCEx            // dcex

} INukeFateHA2025_t;   

class INukeHadroFates2025 {

public:
  //__________________________________________________________________________
  static string AsString(INukeFateHN2025_t fate) {
     switch (fate) {
      case kIHN25FtUndefined : return "** Undefined HN-mode fate **"; break;
      case kIHN25FtCEx       : return "HN-mode / cex";    break;
      case kIHN25FtElas      : return "HN-mode / elas";   break;
      case kIHN25FtInelas    : return "HN-mode / inelas"; break;
      case kIHN25FtAbs       : return "HN-mode / abs";    break;
      case kIHN25FtCmp	   : return "HN-mode / compound"; break;
      case kIHN25FtNoInteraction : return "HN-mode / no interaction"; break;
      default              : break; 
     }
     return "** Undefined HN-mode fate **"; 
  }
  //__________________________________________________________________________
  static string AsString(INukeFateHA2025_t fate) {
     switch (fate) {
      case kIHA25FtUndefined : return "** Undefined HA-mode fate **"; break;
      case kIHA25FtNoInteraction : return "HA-mode / no interaction"; break;
      case kIHA25FtCEx       : return "HA-mode / cex";            break;
	//      case kIHA25FtElas      : return "HA-mode / elas";           break;
      case kIHA25FtInelas    : return "HA-mode / inelas";         break;
      case kIHA25FtAbs       : return "HA-mode / abs";            break;
      case kIHA25FtKo        : return "HA-mode / knock-out";      break;
      case kIHA25FtCmp       : return "HA-mode / compound";       break;
      case kIHA25FtPiProd    : return "HA-mode / pi-production" ; break;
      case kIHA25FtInclPip   : return "HA-mode / pi-prod incl pi+";   break;
      case kIHA25FtInclPim   : return "HA-mode / pi-prod incl pi-";   break;
      case kIHA25FtInclPi0   : return "HA-mode / pi-prod incl pi0";   break;
      case kIHA25FtDCEx      : return "HA-mode / dcex";           break;
      default              : break;
     }
     return "** Undefined HA-mode fate **"; 
  }
  //__________________________________________________________________________
  static string AsSimpleString(INukeFateHA2025_t fate) {
     switch (fate) {
      case kIHA25FtUndefined : return "undefined"; break;
      case kIHA25FtNoInteraction : return "no interaction"; break;
      case kIHA25FtCEx       : return "cex";            break;
	//      case kIHA25FtElas      : return "elas";           break;
      case kIHA25FtInelas    : return "inelas";         break;
      case kIHA25FtAbs       : return "abs";            break;
      case kIHA25FtKo        : return "knock out";      break; 
      case kIHA25FtCmp       : return "compound";       break;
      case kIHA25FtPiProd    : return "pi prod";          break;
      case kIHA25FtDCEx      : return "dcex";             break;
      default              : break;
     }
     return "undefined"; 
  }
  //__________________________________________________________________________

};

}      // genie
#endif // _INTRANUKE_FATES_2025_H_
