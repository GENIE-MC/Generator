//____________________________________________________________________________
/*!

\class    genie::utils::nhl

\brief    Utilities for simulating the decay of Neutral Heavy Leptons 

\author   
          
\created  November 03, 2011

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _NHL_DECAY_UTILS_H_
#define _NHL_DECAY_UTILS_H_

#include <string>

#include "Framework/ParticleData/PDGCodeList.h"
#include "Physics/NeutralHeavyLepton/NHLDecayMode.h"
#include "Physics/NeutralHeavyLepton/NHLProductionMode.h"

// -- for retrieval of parameters from config
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/StringUtils.h"

using std::string;

namespace genie 
{
 namespace utils 
 {
  namespace nhl 
  {
      string       ProdAsString            (genie::NHL::NHLProd_t nhlprod);
      string       AsString                (genie::NHL::NHLDecayMode_t nhldm);
      bool         IsProdKinematicallyAllowed (genie::NHL::NHLProd_t nhldm, double M);
      bool         IsKinematicallyAllowed  (genie::NHL::NHLDecayMode_t nhldm, double Mnhl);
      PDGCodeList  ProductionProductList   (genie::NHL::NHLProd_t nhldm);
      PDGCodeList  DecayProductList        (genie::NHL::NHLDecayMode_t nhldm);

      // for obtaining params, etc, directly from config
      int                 GetCfgInt        (string file_id, string set_name, string par_name);
      std::vector<int>    GetCfgIntVec     (string file_id, string set_name, string par_name);
      double              GetCfgDouble     (string file_id, string set_name, string par_name);
      std::vector<double> GetCfgDoubleVec  (string file_id, string set_name, string par_name);
      bool                GetCfgBool       (string file_id, string set_name, string par_name);
      std::vector<bool>   GetCfgBoolVec    (string file_id, string set_name, string par_name);

   } // nhl 
 } // utils 
} // genie 

#endif // _NHL_DECAY_UTILS_H_
