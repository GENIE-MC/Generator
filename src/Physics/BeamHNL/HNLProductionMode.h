//____________________________________________________________________________
/*!

\class    genie::hnl::HNLProductionMode

\brief    Enumeration of HNL production modes.

\author   John Plows <komninos-john.plows \at physics.ox.ac.uk>

\created  May 06, 2022

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HNL_PRODUCTION_MODE_H_
#define _HNL_PRODUCTION_MODE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {
  namespace hnl {
    
    typedef enum t_HNLProd {
      
      kHNLProdNull          = -1,
      kHNLProdPion2Muon     = 0, // pi --> HNL + mu
      kHNLProdPion2Electron = 1, // pi --> HNL + e
      kHNLProdKaon2Muon     = 2, // K  --> HNL + mu
      kHNLProdKaon2Electron = 3, // K  --> HNL + e
      kHNLProdKaon3Muon     = 4, // K  --> HNL + mu   + pi0
      kHNLProdKaon3Electron = 5, // K  --> HNL + e    + pi0
      kHNLProdNeuk3Muon     = 6, // K0 --> HNL + mu   + pi
      kHNLProdNeuk3Electron = 7, // K0 --> HNL + e    + pi
      kHNLProdMuon3Numu     = 8, // mu --> HNL + numu + e
      kHNLProdMuon3Nue      = 9, // mu --> HNL + nue  + e
      kHNLProdMuon3Nutau    = 10 // mu --> HNL + nutau + e (LFV!)
      
    } HNLProd_t;

  } // namespace hnl
} // namespace genie

#endif // #ifndef _HNL_PRODUCTION_MODE_H_
