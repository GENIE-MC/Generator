//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Julia Tena-Vidal <j.tena-vidal \at liverpool.ac.uk>
 University of Liverpool

*/
//____________________________________________________________________________

#include <cstdlib>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/HadronTransport/CascadeReweight.h"
#include "Physics/NuclearState/NuclearUtils.h"

#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/StringUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//___________________________________________________________________________
CascadeReweight::CascadeReweight()
    : EventRecordVisitorI("genie::CascadeReweight") {}
//___________________________________________________________________________
CascadeReweight::CascadeReweight(string config)
    : EventRecordVisitorI("genie::CascadeReweight", config) {}
//___________________________________________________________________________
CascadeReweight::~CascadeReweight() {}
//___________________________________________________________________________
void CascadeReweight::ProcessEventRecord(GHepRecord *evrec) const {
  if (!evrec) {
    LOG("CascadeReweight", pERROR) << "** Null input!";
    return;
  }
  // Get Associated weight
  double weight = GetEventWeight(*evrec);
  // Set weight
  evrec->SetWeight(weight);

  return;
}
//___________________________________________________________________________
double CascadeReweight::GetEventWeight(const GHepRecord &event) const {

  GHepParticle *p = 0;
  TIter event_iter(&event);
  double total_weight = 1.;
  while ((p = dynamic_cast<GHepParticle *>(event_iter.Next()))) {
    // Look only at stable particles in the nucleus:
    if ( p->Status() != kIStHadronInTheNucleus ) continue ; 
    // Get particle fate
    auto fate_rescatter = p->RescatterCode();
    // Only look at particles that had FSI
    if (fate_rescatter < 0 || fate_rescatter == kIHNFtUndefined)
      continue;
    INukeFateHN_t fate = (INukeFateHN_t)fate_rescatter;

    // Read map weight:
    const auto map_it = fFateWeightsMap.find(fate);
    // Get weight given a pdg code.
    if (map_it != fFateWeightsMap.end()) {
      int pdg_target = p->Pdg();
      const auto weight_it = (map_it->second).find(pdg_target);
      if (weight_it != (map_it->second).end()) {
        total_weight *= weight_it->second;
        continue;
      }
    }
    // If fate is not in the pdg map, use default values:
    const auto def_it = fDefaultMap.find(fate);
    if (def_it != fDefaultMap.end()) {
      total_weight *= def_it->second;
    }
  } // end loop over particles

  return total_weight;
}
//___________________________________________________________________________
void CascadeReweight::Configure(const Registry &config) {
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void CascadeReweight::Configure(string param_set) {
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//___________________________________________________________________________
void CascadeReweight::LoadConfig(void) {
  bool good_config = true;

  // Clean maps
  fDefaultMap.clear();
  fFateWeightsMap.clear();

  // Create vector with list of possible keys (follows the order of the fates
  // enumeration)
  std::map<INukeFateHN_t, string> EINukeFate_map_keys = GetEINukeFateKeysMap();

  for (map<INukeFateHN_t, string>::iterator it_keys =
           EINukeFate_map_keys.begin();
       it_keys != EINukeFate_map_keys.end(); it_keys++) {
    // Find fate specifications
    std::string to_find_def =
        "CascadeReweight-Default-Weight-" + (it_keys->second);

    auto kdef_list = GetConfig().FindKeys(to_find_def.c_str());
    for (auto kiter = kdef_list.begin(); kiter != kdef_list.end(); ++kiter) {
      const RgKey &key = *kiter;
      double weight;
      GetParam(key, weight);
      // Add check weight > 0
      if (weight < 0) {
        LOG("CascadeReweight", pERROR)
            << "The weight assigned to " << to_find_def << " is not positive";
        good_config = false;
        continue;
      }
      fDefaultMap[it_keys->first] = weight;
    }

    // Find Pdg specifications
    std::string to_find_pdg =
        "CascadeReweight-Weight-" + (it_keys->second) + "@Pdg=";
    auto kpdg_list = GetConfig().FindKeys(to_find_pdg.c_str());
    std::map<int, double> WeightMap; // define map that stores <pdg, weight>
    for (auto kiter = kpdg_list.begin(); kiter != kpdg_list.end(); ++kiter) {
      const RgKey &key = *kiter;
      vector<string> kv = genie::utils::str::Split(key, "=");
      assert(kv.size() == 2);
      int pdg_target = stoi(kv[1]);
      if (!PDGLibrary::Instance()->Find(pdg_target)) {
        LOG("CascadeReweight", pERROR)
            << "The target Pdg code associated to " << to_find_pdg
            << " is not valid : " << pdg_target;
        good_config = false;
        continue;
      }
      double weight;
      GetParam(key, weight);
      // Add check weight > 0
      if (weight < 0) {
        LOG("CascadeReweight", pERROR)
            << "The weight assigned to " << to_find_pdg << " is not positive";
        good_config = false;
        continue;
      }
      // Add pdg and weight in map
      WeightMap.insert(std::pair<int, double>(pdg_target, weight));
    }
    // store information in class member
    fFateWeightsMap[it_keys->first] = std::move(WeightMap);
  }

  if (!good_config) {
    LOG("CascadeReweight", pFATAL) << "Configuration has failed.";
    exit(78);
  }
}
//___________________________________________________________________________
