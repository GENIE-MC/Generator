#ifndef _EVTLIB_UTILS_H_
#define _EVTLIB_UTILS_H_

#include <string>

namespace genie{
namespace evtlib{

/// \brief Expand env vars/homedirs/wildcards in \a s
///
/// It is a fatal error if there is not exactly one result of the expansion
void Expand(std::string& s);

}}

#endif
