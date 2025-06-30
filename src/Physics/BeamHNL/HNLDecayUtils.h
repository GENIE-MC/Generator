//____________________________________________________________________________
/*!

\class    genie::utils::hnl

\brief    Utilities for simulating the decay of Heavy Neutral Leptons 

\author   
          
\created  November 03, 2011

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _HNL_DECAY_UTILS_H_
#define _HNL_DECAY_UTILS_H_

#include <string>

#include "Framework/ParticleData/PDGCodeList.h"
#include "Physics/BeamHNL/HNLDecayMode.h"
#include "Physics/BeamHNL/HNLProductionMode.h"

// -- for retrieval of parameters from config
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/StringUtils.h"

using std::string;

namespace genie 
{
 namespace utils 
 {
  namespace hnl 
  {
      string       ProdAsString            (genie::hnl::HNLProd_t hnlprod);
      string       AsString                (genie::hnl::HNLDecayMode_t hnldm);
      bool         IsProdKinematicallyAllowed (genie::hnl::HNLProd_t hnlprod);
      bool         IsKinematicallyAllowed  (genie::hnl::HNLDecayMode_t hnldm, double Mhnl);
      PDGCodeList  ProductionProductList   (genie::hnl::HNLProd_t hnldm);
      PDGCodeList  DecayProductList        (genie::hnl::HNLDecayMode_t hnldm);

      // for obtaining params, etc, directly from config
      int                 GetCfgInt        (string file_id, string set_name, string par_name);
      std::vector<int>    GetCfgIntVec     (string file_id, string set_name, string par_name);
      double              GetCfgDouble     (string file_id, string set_name, string par_name);
      std::vector<double> GetCfgDoubleVec  (string file_id, string set_name, string par_name);
      bool                GetCfgBool       (string file_id, string set_name, string par_name);
      std::vector<bool>   GetCfgBoolVec    (string file_id, string set_name, string par_name);
      std::string         GetCfgString     (string file_id, string set_name, string par_name);

   } // hnl 
 } // utils 
} // genie 

#endif // _HNL_DECAY_UTILS_H_
