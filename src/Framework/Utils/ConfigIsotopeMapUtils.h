//____________________________________________________________________________
/*!

\namespace  genie::utils::config

\brief      Simple functions for loading and reading nucleus dependent keys
            from config files.

\author     Brian Coopersmith, University of Rochester

\created    October 23, 2014

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef _CONFIG_UTILS_H
#define _CONFIG_UTILS_H

#include <map>
#include "Framework/Interaction/Target.h"

namespace genie {
class Registry;

namespace utils {
namespace config {
bool GetValueFromNuclearMaps(
    const Target& target, const std::map<int, double>& nuc_to_val,
    const std::map<std::pair<int, int>, double>& nucA_range_to_val,
    double* val);
void LoadAllNucARangesForKey(const char* key_name, const char* log_tool_name,
                             Registry* config,
                             std::map<std::pair<int, int>, double>* nuc_rangeA_to_val);
void LoadAllIsotopesForKey(const char* key_name, const char* log_tool_name,
                           Registry* config, std::map<int, double>* nuc_to_val);

bool GetDoubleKeyPDG(const char* valName, const int pdgc,
                     Registry* config, double* val);
bool GetDoubleKeyRangeNucA(const char* valName, const int lowA,
                           const int highA, Registry* config, double* val);

}  // namespace config
}  // namespace utils
}  // namespace genie

#endif  // _CONFIG_UTILS_H
