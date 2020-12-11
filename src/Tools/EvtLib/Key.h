#ifndef _EVTLIB_KEY_H_
#define _EVTLIB_KEY_H_

#include <ostream>

namespace genie{
namespace evtlib{

struct Key
{
  Key(int _nucl_pdg, int _nu_pdg, bool _iscc)
    : nucl_pdg(_nucl_pdg), nu_pdg(_nu_pdg), iscc(_iscc) {}

  bool operator<(const Key& k) const
  {
    return (std::make_tuple(  nucl_pdg,   nu_pdg,   iscc) <
            std::make_tuple(k.nucl_pdg, k.nu_pdg, k.iscc));
  }

  friend std::ostream& operator<<(std::ostream& os, const Key& k)
  {
    os << k.nu_pdg << " on " << k.nucl_pdg << " " << " via " << (k.iscc ? "CC" : "NC");
    return os;
  }

  int nucl_pdg;
  int nu_pdg;
  bool iscc;
};

}} // namespaces

#endif
