#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCTreeMultiHeader.h"
#include "Framework/Messenger/Messenger.h"

#include <algorithm>
#include <set>

using namespace genie;

std::string NtpMCTreeMultiHeader::getFileName(size_t i) const
{
    return checkSize(i) ? inputfiles[i] : "";
}

std::pair<Long64_t, Long64_t> NtpMCTreeMultiHeader::getIndices(size_t i) const
{
    return checkSize(i) ? indices[i] : std::make_pair(Long64_t{-1}, Long64_t{-1});
}

std::set<int> NtpMCTreeMultiHeader::getNeutrinos(size_t i) const
{
    return checkSize(i) ? neutrinos[i] : std::set<int>{};
}

const NtpMCTreeHeader NtpMCTreeMultiHeader::getHead(size_t i) const
{
    return checkSize(i) ? heads[i] : NtpMCTreeHeader{};
}

bool NtpMCTreeMultiHeader::checkSize(size_t i) const
{
    const bool valid = i < size();

    if (!valid) {
        LOG("NtpMCTreeMultiHeader", pWARN)
            << "Requesting nonexistent Head, ignoring";
    }

    return valid;
}

size_t NtpMCTreeMultiHeader::size() const
{
    return std::min({
        inputfiles.size(),
        indices.size(),
        neutrinos.size(),
        heads.size()
    });
}

void NtpMCTreeMultiHeader::insertHead(
    NtpMCTreeHeader head,
    std::string filename,
    std::set<int> neutrino,
    std::pair<Long64_t, Long64_t> indexrange
    )
{
    inputfiles.push_back(std::move(filename));
    heads.push_back(std::move(head));
    indices.push_back(indexrange);
    neutrinos.push_back(std::move(neutrino));
}