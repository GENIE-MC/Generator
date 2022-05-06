//____________________________________________________________________________
/*!

\class    genie::NHL::NHLProductionMode

\brief    Enumeration of NHL production modes.

\author   John Plows <komninos-john.plows \at physics.ox.ac.uk>

\created  May 06, 2022

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NHL_PRODUCTION_MODE_H_
#define _NHL_PRODUCTION_MODE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {
  namespace NHL {
    
    typedef enum t_NHLProd {
      
      kNHLProdPion2Muon     = 0, // pi --> NHL + mu
      kNHLProdPion2Electron = 1, // pi --> NHL + e
      kNHLProdKaon2Muon     = 2, // K  --> NHL + mu
      kNHLProdKaon2Electron = 3, // K  --> NHL + e
      kNHLProdKaon3Muon     = 4, // K  --> NHL + mu   + pi0
      kNHLProdKaon3Electron = 5, // K  --> NHL + e    + pi0
      kNHLProdNeuk3Muon     = 6, // K0 --> NHL + mu   + pi
      kNHLProdNeuk3Electron = 7, // K0 --> NHL + e    + pi
      kNHLProdMuon3Numu     = 8, // mu --> NHL + numu + e
      kNHLProdMuon3Nue      = 9, // mu --> NHL + nue  + e
      kNHLProdMuon3Nutau    = 10 // mu --> NHL + nutau + e (LFV!)
      
    } NHLProd_t;

  } // namespace NHL
} // namespace genie

#endif // #ifndef _NHL_PRODUCTION_MODE_H_
