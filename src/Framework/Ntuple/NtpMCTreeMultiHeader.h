#ifndef MULTIHEAD_H
#define MULTIHEAD_H
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include <vector>
#include <string>
#include <utility>


namespace genie {
class NtpMCTreeMultiHeader : public TObject {
public:
NtpMCTreeMultiHeader(){};
~NtpMCTreeMultiHeader()=default;
  std::vector<std::string> inputfiles;

  // [first event, last event] in the merged output corresponding
  // to each input file.
  std::vector<std::pair<Long64_t, Long64_t>> indices;
  // generated neutrinos each run.
  std::vector<std::vector<int>> neutrinos;
  // Store the headers of all input files.
  std::vector<NtpMCTreeHeader> heads;
  ClassDef(NtpMCTreeMultiHeader, 1)
};
};

#endif