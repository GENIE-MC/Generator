////////////////////////////////////////////////////////////////////////
// \author Christopher Backhouse -- c.backhouse@ucl.ac.uk
////////////////////////////////////////////////////////////////////////

#include "Tools/VMC/RecordList.h"

#include "Framework/Messenger/Messenger.h"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>

namespace genie{
namespace vmc{
  //---------------------------------------------------------------------------
  Record::Record() : E(0)/*, weight(0)*/
  {
  }

  //---------------------------------------------------------------------------
    Record::Record(float _E, /*float _w,*/ int _prod_id, const std::vector<Particle>& _ps)
      : E(_E), /*weight(_w),*/ prod_id(_prod_id), parts(_ps)
  {
  }

  //---------------------------------------------------------------------------
  bool Record::operator<(const Record& rhs) const
  {
    return E < rhs.E;
  }

  //---------------------------------------------------------------------------
  RecordLoader::RecordLoader(const char* fname)
  {
    fFile = new TFile(fname);
    if(fFile->IsZombie()) exit(1);

    fTree = (TTree*)fFile->Get("tr");
    if(!fTree){
      LOG("ELI", pFATAL) << "'tr' not found in " << fname;
      exit(1);
    }

    fTree->SetBranchAddress("Enu", &Enu);
//    fTree->SetBranchAddress("weight", &weight);
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
    delete fFile;
  }

  //---------------------------------------------------------------------------
  long RecordLoader::NRecords() const
  {
    return fTree->GetEntries();
  }

  //---------------------------------------------------------------------------
  Record RecordLoader::GetRecord(int i) const
  {
    fTree->GetEntry(i);

    if(nparts > 1024){
      LOG("ELI", pFATAL) << "Too many particles " << nparts;
      exit(1);
    }

    std::vector<Particle> parts(nparts);
    for(int j = 0; j < nparts; ++j){
      parts[j].pdg = pdgs[j];
      parts[j].E = Es[j];
      parts[j].px = px[j];
      parts[j].py = py[j];
      parts[j].pz = pz[j];
    } // end for j

    return Record(Enu, /*weight,*/ prod_id, parts);
  }

  //---------------------------------------------------------------------------
  SimpleRecordList::SimpleRecordList(const char* fname)
  {
    std::cout << "Loading " << fname;
    RecordLoader loader(fname);

    const int N = loader.NRecords();
    fRecs.reserve(N);
    for(int i = 0; i < N; ++i){
      if(i%(N/8) == 0) std::cout << "." << std::flush;

      fRecs.push_back(loader.GetRecord(i));
    } // end for i
    std::cout << std::endl;

    std::sort(fRecs.begin(), fRecs.end());

//    for(Record& r: fRecs) r.weight *= fRecs.size();
  }

  //---------------------------------------------------------------------------
  const Record* SimpleRecordList::GetRecord(float E) const
  {
    auto it = std::lower_bound(fRecs.begin(), fRecs.end(), Record(E, /*0,*/ 0, {}));
    if(it == fRecs.end()) return 0;
    return &(*it);
  }

  //---------------------------------------------------------------------------
  OnDemandRecordList::OnDemandRecordList(const char* fname)
    : fFname(fname), fLoader(fname)
  {
  }

  //---------------------------------------------------------------------------
  void OnDemandRecordList::LoadIndex() const
  {
    std::cout << "Loading index to " << fFname;

    TFile f(fFname.c_str());
    if(f.IsZombie()) exit(1);

    TTree* tr = (TTree*)f.Get("tr");
    if(!tr){
      LOG("ELI", pFATAL) << "'tr' not found in " << fFname;
      exit(1);
    }

    float Enu;
    tr->SetBranchAddress("Enu", &Enu);

    const int N = tr->GetEntries();
    fEnergies.reserve(N);

    for(int i = 0; i < N; ++i){
      if(i%(N/8) == 0) std::cout << "." << std::flush;
      tr->GetEntry(i);

      fEnergies.emplace_back(Enu, i);
    } // end for i
    std::cout << std::endl;

    std::sort(fEnergies.begin(), fEnergies.end());
  }

  //---------------------------------------------------------------------------
  const Record* OnDemandRecordList::GetRecord(float E) const
  {
    if(fEnergies.empty()) LoadIndex();

    auto it = std::lower_bound(fEnergies.begin(), fEnergies.end(),
                               std::make_pair(E, 0));
    if(it == fEnergies.end()) return 0;

    fRecord = fLoader.GetRecord(it->second);

//    fRecord.weight *= fEnergies.size();

    return &fRecord;
  }
}} // namespaces
