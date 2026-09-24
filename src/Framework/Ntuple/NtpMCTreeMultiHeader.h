#ifndef MULTIHEAD_H
#define MULTIHEAD_H

#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "TFolder.h"
#include "TObject.h"

#include <set>
#include <string>
#include <utility>
#include <vector>

namespace genie {

class NtpMCTreeMultiHeader : public TObject {
  public:
    NtpMCTreeMultiHeader() = default;
    ~NtpMCTreeMultiHeader() = default;

    size_t size() const;

    std::string getFileName(size_t i) const;

    std::pair<Long64_t, Long64_t> getIndices(size_t i) const;

    std::set<int> getNeutrinos(size_t i) const;

    const NtpMCTreeHeader getHead(size_t i) const;
    void PrintToStream(ostream & stream) const;

    void insertHead(
        NtpMCTreeHeader head,
        std::string filename,
        std::set<int> neutrino,
        std::pair<Long64_t, Long64_t> indexrange
    );

private:
    bool checkSize(size_t i) const;

    std::vector<std::string> inputfiles;

    // [first event, last event] in the merged output corresponding
    // to each input file.
    std::vector<std::pair<Long64_t, Long64_t>> indices;

    // Generated neutrinos each run.
    std::vector<std::set<int>> neutrinos;

    // Store the headers of all input files.
    std::vector<NtpMCTreeHeader> heads;


    ClassDef(NtpMCTreeMultiHeader, 1)
};
ostream & operator << (ostream & stream, const NtpMCTreeMultiHeader & hdr);

} // namespace genie

#endif
