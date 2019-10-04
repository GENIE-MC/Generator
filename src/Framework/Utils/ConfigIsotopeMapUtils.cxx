//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Brian Coopersmith, University of Rochester

 For documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/Utils/ConfigIsotopeMapUtils.h"

#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"

using namespace std;

namespace genie {
namespace utils {
namespace config {

//____________________________________________________________________________
// Given a map of nucleus PDGs to values and ranges of nucleus As to values,
// return the value for the supplied target. Individual nucleus maps take 
// precedence over the range maps.
//____________________________________________________________________________
bool GetValueFromNuclearMaps(
    const Target& target, const map<int, double>& nuc_to_val,
    const map<pair<int, int>, double>& nucA_range_to_val,
    double* val) {
  const int pdgc = pdg::IonPdgCode(target.A(), target.Z());
  map<int, double>::const_iterator nuc_it = nuc_to_val.find(pdgc);
  if(nuc_it != nuc_to_val.end()) {
    *val = nuc_it->second;
    return true;
  }
  map<pair<int, int>, double>::const_iterator range_it =
      nucA_range_to_val.begin();
  for(; range_it != nucA_range_to_val.end(); ++range_it) {
    if (target.A() >= range_it->first.first &&
        target.A() <= range_it->first.second) {
      *val = range_it->second;
      return true;
    }
  }
  return false;
}
//____________________________________________________________________________
// Read in from the config file all listed NucA range parameters for a given
// key.  Valid for As up to 419.  Stores them in the map from
// pair(lowA, highA) to value.
//____________________________________________________________________________
void LoadAllNucARangesForKey(const char* key_name, const char* log_tool_name,
                             Registry* config,
                             map<pair<int, int>, double>* nuc_rangeA_to_val) {
  for (int lowA = 1; lowA < 3 * 140; lowA++) {
    for (int highA = lowA; highA < 3 * 140; highA++) {
      double val;
      if (GetDoubleKeyRangeNucA(key_name, lowA, highA, config, &val)) {
        LOG(log_tool_name, pINFO) << "For "<< lowA - 1 <<" < A < " <<
            highA + 1 << " -> using " << key_name << " = " << val;
        (*nuc_rangeA_to_val)[pair<int, int>(lowA, highA)] = val;
      }
    }
  }
}
//____________________________________________________________________________
// Read in from the config file all listed NucZ range parameters for a given
// key. Valid for Zs up to 139 and As up to 3*Z.  Stores them in the map from
// PDG code to value.
//____________________________________________________________________________
void LoadAllIsotopesForKey(const char* key_name, const char* log_tool_name,
                           Registry* config, map<int, double>* nuc_to_val) {
  for (int Z = 1; Z < 140; Z++) {
    for (int A = Z; A < 3 * Z; A++) {
      const int pdgc = pdg::IonPdgCode(A, Z);
      double val;
      if(GetDoubleKeyPDG(key_name, pdgc, config, &val)) {
        LOG(log_tool_name, pINFO) << "Nucleus: " << pdgc <<
            " -> using " << key_name << " = " << val;
        (*nuc_to_val)[pdgc] = val;
      }
    }
  }
}
//____________________________________________________________________________
// Check if the key <valName>@Pdg=<pdgc> exists in config.  If so, load that
// into val, and return true.  Otherwise return false.
//____________________________________________________________________________
bool GetDoubleKeyPDG(const char* valName, const int pdgc,
                     Registry* config, double* val)
{
  ostringstream s;
  s<<valName<<"@Pdg="<<pdgc;
  RgKey key = s.str();
  if(!config->Exists(key)) {
    return false;
  }
  *val = config->GetDoubleDef(key,0);
  return true;
}
//____________________________________________________________________________
// Check if the key <valName>@LowA=<lowA>;HighA=<highA> exists in config. If
// so load that into val and return true.  Otherwise return false.
//____________________________________________________________________________
bool GetDoubleKeyRangeNucA(const char* valName, const int lowA, 
                           const int highA, Registry* config, double* val)
{
  ostringstream s;
  s<<valName<<"@LowA="<<lowA<<";HighA="<<highA;
  RgKey key = s.str();
  if(!config->Exists(key)) {
    return false;
  }
  *val = config->GetDoubleDef(key,0);
  return true;
}

}  // namespace config
}  // namespace utils
}  // namespace genie
