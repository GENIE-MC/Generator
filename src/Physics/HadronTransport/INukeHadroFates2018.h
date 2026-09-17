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

#ifndef _INTRANUKE_FATES_2018_H_
#define _INTRANUKE_FATES_2018_H_

#include <string>

using std::string;

namespace genie {

// Fates in INTRANUKE's HN mode
//
typedef enum EINukeFateHN2018_t {

   kIHN18FtUndefined = 0, 
   kIHN18FtNoInteraction,
   kIHN18FtCEx,       // cex
   kIHN18FtElas,      // elas
   kIHN18FtInelas,    // inelas
   kIHN18FtAbs,       // abs 
   kIHN18FtCmp         //cmp

} INukeFateHN2018_t;   

// Fates in INTRANUKE's HA mode
//
typedef enum EINukeFateHA2018_t {

   kIHA18FtUndefined = 0,
   kIHA18FtNoInteraction,  // no interaction 
   kIHA18FtCEx,            // cex
   //   kIHA18FtElas,           // elas
   kIHA18FtInelas,         // inelas
   kIHA18FtAbs,            // abs
   kIHA18FtKo, 	         // knock out
   kIHA18FtCmp,            // compound nucleus
   kIHA18FtPiProd,         // pi production
   kIHA18FtInclPip,        // pi production : inclusive pi+
   kIHA18FtInclPim,        // pi production : inclusive pi-
   kIHA18FtInclPi0,        // pi production : inclusive pi0 
   kIHA18FtDCEx            // dcex

} INukeFateHA2018_t;   

class INukeHadroFates2018 {

public:
  //__________________________________________________________________________
  static string AsString(INukeFateHN2018_t fate) {
     switch (fate) {
      case kIHN18FtUndefined : return "** Undefined HN-mode fate **"; break;
      case kIHN18FtCEx       : return "HN-mode / cex";    break;
      case kIHN18FtElas      : return "HN-mode / elas";   break;
      case kIHN18FtInelas    : return "HN-mode / inelas"; break;
      case kIHN18FtAbs       : return "HN-mode / abs";    break;
      case kIHN18FtCmp	   : return "HN-mode / compound"; break;
      case kIHN18FtNoInteraction : return "HN-mode / no interaction"; break;
      default              : break; 
     }
     return "** Undefined HN-mode fate **"; 
  }
  //__________________________________________________________________________
  static string AsString(INukeFateHA2018_t fate) {
     switch (fate) {
      case kIHA18FtUndefined : return "** Undefined HA-mode fate **"; break;
      case kIHA18FtNoInteraction : return "HA-mode / no interaction"; break;
      case kIHA18FtCEx       : return "HA-mode / cex";            break;
	//      case kIHA18FtElas      : return "HA-mode / elas";           break;
      case kIHA18FtInelas    : return "HA-mode / inelas";         break;
      case kIHA18FtAbs       : return "HA-mode / abs";            break;
      case kIHA18FtKo        : return "HA-mode / knock-out";      break;
      case kIHA18FtCmp       : return "HA-mode / compound";       break;
      case kIHA18FtPiProd    : return "HA-mode / pi-production" ; break;
      case kIHA18FtInclPip   : return "HA-mode / pi-prod incl pi+";   break;
      case kIHA18FtInclPim   : return "HA-mode / pi-prod incl pi-";   break;
      case kIHA18FtInclPi0   : return "HA-mode / pi-prod incl pi0";   break;
      case kIHA18FtDCEx      : return "HA-mode / dcex";           break;
      default              : break;
     }
     return "** Undefined HA-mode fate **"; 
  }
  //__________________________________________________________________________
  static string AsSimpleString(INukeFateHA2018_t fate) {
     switch (fate) {
      case kIHA18FtUndefined : return "undefined"; break;
      case kIHA18FtNoInteraction : return "no interaction"; break;
      case kIHA18FtCEx       : return "cex";            break;
	//      case kIHA18FtElas      : return "elas";           break;
      case kIHA18FtInelas    : return "inelas";         break;
      case kIHA18FtAbs       : return "abs";            break;
      case kIHA18FtKo        : return "knock out";      break; 
      case kIHA18FtCmp       : return "compound";       break;
      case kIHA18FtPiProd    : return "pi prod";          break;
      case kIHA18FtDCEx      : return "dcex";             break;
      default              : break;
     }
     return "undefined"; 
  }
  //__________________________________________________________________________

};

}      // genie
#endif // _INTRANUKE_FATES_2018_H_
