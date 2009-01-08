//____________________________________________________________________________
/*!

\class    genie::INukeHadroFates

\brief    An enumeration of possible hadron "fates" taken into account by the
          INTRANUKE hadron transport MC.

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>, Rutherford Lab.

\created  November 1, 2005

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE 
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
   kIHNFtAbs        // abs 

} INukeFateHN_t;   

// Fates in INTRANUKE's HA mode
//
typedef enum EINukeFateHA_t {

   kIHAFtUndefined = 0, 
   kIHAFtCEx,       // cex
   kIHAFtElas,      // elas
   kIHAFtInelas,    // inelas
   kIHAFtAbsNP,     // abs np
   kIHAFtAbsPP,     // abs pp
   kIHAFtAbsNPP,    // abs npp
   kIHAFtAbsNNP,    // abs nnp
   kIHAFtAbs2N2P,   // abs 2n2p
   kIHAFtAbs2N3P,   // abs 2n3p
   kIHAFtNPip,      // pi production : n pi+
   kIHAFtNPipPi0    // pi production : n pi+ pi0

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
      default              : break; 
     }
     return "** Undefined HN-mode fate **"; 
  }
  //__________________________________________________________________________
  static string AsString(INukeFateHA_t fate) {
     switch (fate) {
      case kIHAFtUndefined : return "** Undefined HA-mode fate **"; break;
      case kIHAFtCEx       : return "HA-mode / cex";          break;
      case kIHAFtElas      : return "HA-mode / elas";         break;
      case kIHAFtInelas    : return "HA-mode / inelas";       break;
      case kIHAFtAbsNP     : return "HA-mode / abs pn";       break;
      case kIHAFtAbsPP     : return "HA-mode / abs pp";       break;
      case kIHAFtAbsNPP    : return "HA-mode / abs npp";      break;
      case kIHAFtAbsNNP    : return "HA-mode / abs nnp";      break;
      case kIHAFtAbs2N2P   : return "HA-mode / abs 2n2p";     break;
      case kIHAFtAbs2N3P   : return "HA-mode / abs 2n3p";     break;
      case kIHAFtNPip      : return "HA-mode / abs npi+";     break;
      case kIHAFtNPipPi0   : return "HA-mode / abs npi+pi0";  break;
      default              : break;
     }
     return "** Undefined HA-mode fate **"; 
  }
  //__________________________________________________________________________
};

}      // genie
#endif // _INTRANUKE_FATES_H_
