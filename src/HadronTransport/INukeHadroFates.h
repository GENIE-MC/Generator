//____________________________________________________________________________
/*!

\class    genie::INukeHadroFates

\brief    An enumeration of possible hadron "fates" taken into account by the
          INTRANUKE hadron transport MC.

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
          Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, Rutherford Lab.

\created  November 1, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE. 
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_FATES_H_
#define _INTRANUKE_FATES_H_

#include <string>

using std::string;

namespace genie {

// Fates in INTRANUKE's HN mode
//
typedef enum EINukeFateHN_t {

   kIHNFtUndefined = 0, 
   kIHNFtCEx,       // cex
   kIHNFtElas,      // elas
   kIHNFtInelas,    // inelas
   kIHNFtAbsPN      // abs 

} INukeFateHN_t;   

// Fates in INTRANUKE's HA mode
//
typedef enum EINukeFateHA_t {

   kIHAFtUndefined = 0, 
   kIHAFtCEx,       // cex
   kIHAFtElas,      // elas
   kIHAFtInelas,    // inelas
   kIHAFtAbsPN,     // abs pn
   kIHAFtAbsPP,     // abs pp
   kIHAFtAbsNPP,    // abs npp
   kIHAFtAbsNNP,    // abs nnp
   kIHAFtAbs4N4P,   // abs 4n4p (remainder)
   kIHAFtPiProd,    // piprod
   kIHAFt10         // fate 10 ??

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
      case kIHNFtAbsPN     : return "HN-mode / abs";    break;
      default              : break; 
     }
     return "** Undefined HN-mode fate **"; 
  }
  //__________________________________________________________________________
  static string AsString(INukeFateHA_t fate) {
     switch (fate) {
      case kIHAFtUndefined : return "** Undefined HA-mode fate **"; break;
      case kIHAFtCEx       : return "HA-mode / cex";      break;
      case kIHAFtElas      : return "HA-mode / elas";     break;
      case kIHAFtInelas    : return "HA-mode / inelas";   break;
      case kIHAFtAbsPN     : return "HA-mode / abs pn";   break;
      case kIHAFtAbsPP     : return "HA-mode / abs pp";   break;
      case kIHAFtAbsNPP    : return "HA-mode / abs npp";  break;
      case kIHAFtAbsNNP    : return "HA-mode / ans nnp";  break;
      case kIHAFtAbs4N4P   : return "HA-mode / abs 4n4p"; break;
      case kIHAFtPiProd    : return "HA-mode / piprod";   break;
      case kIHAFt10        : return "HA-mode / fate 10?"; break;
      default              : break;
     }
     return "** Undefined HA-mode fate **"; 
  }
  //__________________________________________________________________________
};

}      // genie
#endif // _INTRANUKE_FATES_H_
