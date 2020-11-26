////////////////////////////////////////////////////////////////////////
// \author Chris Backhouse -- c.backhouse@ucl.ac.uk
////////////////////////////////////////////////////////////////////////

#include "Tools/EvtLib/EvtLibRecordList.h"

#include "Framework/Messenger/Messenger.h"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>

namespace genie{
namespace evtlib{
  //---------------------------------------------------------------------------
  EvtLibRecord::EvtLibRecord() : E(0)
  {
  }

  //---------------------------------------------------------------------------
  EvtLibRecord::EvtLibRecord(float _E, int _prod_id,
                             const std::vector<EvtLibParticle>& _ps)
    : E(_E), prod_id(_prod_id), parts(_ps)
  {
  }

  //---------------------------------------------------------------------------
  bool EvtLibRecord::operator<(const EvtLibRecord& rhs) const
  {
    return E < rhs.E;
  }

  //---------------------------------------------------------------------------
  RecordLoader::RecordLoader(TTree* tree)
    : fTree(tree)
  {
    fTree->SetBranchAddress("Enu", &Enu);
    fTree->SetBranchAddress("prod_id", &prod_id);

    fTree->SetBranchAddress("nparts", &nparts);

    fTree->SetBranchAddress("pdg", &pdgs);
    fTree->SetBranchAddress("E", &Es);
    fTree->SetBranchAddress("px", &px);
    fTree->SetBranchAddress("py", &py);
    fTree->SetBranchAddress("pz", &pz);
  }

  //---------------------------------------------------------------------------
  RecordLoader::~RecordLoader()
  {
  }

  //---------------------------------------------------------------------------
  long RecordLoader::NRecords() const
  {
    return fTree->GetEntries();
  }

  //---------------------------------------------------------------------------
  EvtLibRecord RecordLoader::GetRecord(int i) const
  {
    fTree->GetEntry(i);

    // The event library ROOT format is as minimalistic as possible. We use
    // variable-sized arrays for the particle list. Because the address of
    // `parts` had to be provided to SetBranchAddress up-front we had to pick
    // a fixed size. If the list was actually longer we're in trouble and
    // should bail out. This seems extremely unlikely to be a problem in
    // practice, and we can always increase the constant to something larger.
    if(nparts > kEvtLibMaxParts){
      LOG("ELI", pFATAL) << "Too many particles " << nparts << "(limit " << kEvtLibMaxParts << ")";
      exit(1);
    }

    std::vector<EvtLibParticle> parts(nparts);
    for(int j = 0; j < nparts; ++j){
      parts[j].pdg = pdgs[j];
      parts[j].E = Es[j];
      parts[j].px = px[j];
      parts[j].py = py[j];
      parts[j].pz = pz[j];
    } // end for j

    return EvtLibRecord(Enu, prod_id, parts);
  }

  //---------------------------------------------------------------------------
  SimpleRecordList::SimpleRecordList(TTree* tree, const std::string& prettyName)
  {
    std::cout << "Loading " << prettyName;
    RecordLoader loader(tree);

    const int N = loader.NRecords();
    fRecs.reserve(N);
    for(int i = 0; i < N; ++i){
      if(i%(N/8) == 0) std::cout << "." << std::flush;

      fRecs.push_back(loader.GetRecord(i));
    } // end for i
    std::cout << std::endl;

    std::sort(fRecs.begin(), fRecs.end());
  }

  //---------------------------------------------------------------------------
  const EvtLibRecord* SimpleRecordList::GetRecord(float E) const
  {
    auto it = std::lower_bound(fRecs.begin(), fRecs.end(),
                               EvtLibRecord(E, 0, {}));
    if(it == fRecs.end()) return 0;
    return &(*it);
  }

  //---------------------------------------------------------------------------
  OnDemandRecordList::OnDemandRecordList(TTree* tree, const std::string& prettyName)
    : fTree(tree), fPrettyName(prettyName), fLoader(0)
  {
  }

  //---------------------------------------------------------------------------
  void OnDemandRecordList::LoadIndex() const
  {
    std::cout << "Loading index to " << fPrettyName;

    float Enu;
    fTree->SetBranchAddress("Enu", &Enu);

    const int N = fTree->GetEntries();
    fEnergies.reserve(N);

    for(int i = 0; i < N; ++i){
      if(i%(N/8) == 0) std::cout << "." << std::flush;
      fTree->GetEntry(i);

      fEnergies.emplace_back(Enu, i);
    } // end for i
    std::cout << std::endl;

    std::sort(fEnergies.begin(), fEnergies.end());
  }

  //---------------------------------------------------------------------------
  const EvtLibRecord* OnDemandRecordList::GetRecord(float E) const
  {
    if(fEnergies.empty()){
      LoadIndex();
      // The loader must be created after the indices are loaded, because they
      // both want to set branch addresses in the TTree, and we ultimately need
      // the loader's to be active.
      fLoader = new RecordLoader(fTree);
    }

    auto it = std::lower_bound(fEnergies.begin(), fEnergies.end(),
                               std::make_pair(E, 0));
    if(it == fEnergies.end()) return 0;

    fRecord = fLoader->GetRecord(it->second);

    return &fRecord;
  }
}} // namespaces
