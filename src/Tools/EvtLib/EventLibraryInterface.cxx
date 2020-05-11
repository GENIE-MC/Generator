//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org


*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Interaction/Interaction.h"
#include "Tools/EvtLib/EventLibraryInterface.h"
#include "Tools/EvtLib/EvtLibRecordList.h"
#include "Tools/EvtLib/Utils.h"

#include "TFile.h"

using namespace genie;
using namespace genie::evtlib;

//____________________________________________________________________________
EventLibraryInterface::EventLibraryInterface() :
  EventRecordVisitorI("genie::evtlib::EventLibraryInterface"),
  fRecordFile(0)
{

}

//____________________________________________________________________________
EventLibraryInterface::EventLibraryInterface(string config) :
  EventRecordVisitorI("genie::evtlib::EventLibraryInterface", config),
  fRecordFile(0)
{

}

//____________________________________________________________________________
EventLibraryInterface::~EventLibraryInterface()
{
  Cleanup();
}

//____________________________________________________________________________
void EventLibraryInterface::ProcessEventRecord(GHepRecord * event) const
{
// Get event summary constructed by GENIE
//
  Interaction* interaction = event->Summary();
  const InitialState & init_state = interaction->InitState();

  const EvtLibRecord* rec = GetRecord(interaction);
  if(!rec) return; // Reason has already been printed

  std::unique_ptr<TLorentzVector> probe_p4(init_state.GetProbeP4(kRfLab));

  unsigned int nLep = 0;
  for(const EvtLibParticle& part: rec->parts) if(pdg::IsLepton(part.pdg)) ++nLep;

  const int firstLep = (nLep == 0) ? -1 : 2;
  const int lastLep  = (nLep == 0) ? -1 : 1+nLep;

  const int firstHad = (rec->parts.size() == nLep) ? -1 : 2+nLep;
  const int lastHad  = (rec->parts.size() == nLep) ? -1 : 1+rec->parts.size();

  // Neutrino is a parent to the lepton(s)
  event->AddParticle(init_state.ProbePdg(),
                     kIStInitialState,
                     -1, -1, firstLep, lastLep,
                     *probe_p4,
                     TLorentzVector(0, 0, 0, 0));

  const int tgt_A    = init_state.Tgt().A();
  const int tgt_Z    = init_state.Tgt().Z();
  const int tgt_pdgc = pdg::IonPdgCode(tgt_A, tgt_Z);
  // REVIEW: it was like this in Costas' example code. Can it just be
  // init_state.TgtPdg()?

  LOG("ELI", pINFO)
    << "Adding nucleus [A = " << tgt_A << ", Z = " << tgt_Z
    << ", pdg = " << tgt_pdgc << "]";

  // Nucleus is a parent to everything else
  event->AddParticle(tgt_pdgc,
                     kIStInitialState,
                     -1, -1, firstHad, lastHad,
                     0, 0, 0, PDGLibrary::Instance()->Find(tgt_pdgc)->Mass(),
                     0, 0, 0, 0);

  const std::vector<TVector3> basis = Basis(probe_p4->Vect());

  for(bool lep: {true, false}){ // make sure lepton in first
    for(const EvtLibParticle& part: rec->parts){
      if(lep != pdg::IsLepton(part.pdg)) continue;

      // Fix up the outgoing lepton for NC events (in the library it's always
      // nu_mu...)
      int pdg = part.pdg;
      if(interaction->ProcInfo().IsWeakNC()) pdg = init_state.ProbePdg();

      event->AddParticle(pdg,
                         kIStStableFinalState,
                         (lep ? 0 : 1), -1, -1, -1, // child of the neutrino or nucleus
                         TLorentzVector(part.px*basis[0] +
                                        part.py*basis[1] +
                                        part.pz*basis[2],
                                        part.E),
                         TLorentzVector(0, 0, 0, 0));
    }
  }
}

//____________________________________________________________________________
const EvtLibRecord* EventLibraryInterface::
GetRecord(const Interaction* interaction) const
{
  const InitialState& init_state = interaction->InitState();

  const double probe_E = init_state.ProbeE(kRfLab);

  if(!init_state.Tgt().IsNucleus()){
    LOG("ELI", pINFO) << "Skippping non-nuclear target " << init_state;
    return 0;
  }

  const int tgt_pdgc = init_state.TgtPdg();

  const ProcessInfo& proc = interaction->ProcInfo();

  if(!proc.IsWeakCC() && !proc.IsWeakNC()){
    LOG("ELI", pINFO) << "Skipping unknown process " << proc;
    return 0;
  }

  int probe_pdgc = init_state.ProbePdg();

  // Use nu_mu for NC as a convention internal to this code to index into the
  // records map.
  if(proc.IsWeakNC()){
    if(probe_pdgc > 0) probe_pdgc = +14; else probe_pdgc = -14;
  }

  const Key key(tgt_pdgc, probe_pdgc, proc.IsWeakCC());

  const auto rec_it = fRecords.find(key);

  if(rec_it == fRecords.end()){
    LOG("ELI", pINFO) << "Skippping " << key << " -- not found in library";
    return 0;
  }

  const EvtLibRecord* rec = rec_it->second->GetRecord(probe_E);

  if(!rec){
    LOG("ELI", pINFO) << "Skippping " << key << " at " << probe_E << " GeV -- not found in library";
    return 0;
  }

  return rec;
}

//____________________________________________________________________________
std::vector<TVector3> EventLibraryInterface::Basis(const TVector3& z) const
{
  TVector3 up(0, 1, 0);
  if(up.Dot(z) == 0) up = TVector3(1, 0, 0);

  const TVector3 x = up.Cross(z).Unit(); // Perpendicular to neutrino and up
  const TVector3 y = x.Cross(z).Unit();  // Defines the third axis

  const double a = RandomGen::Instance()->RndEvg().Uniform(0, 2*M_PI);

  const TVector3 xp =  cos(a) * x + sin(a) * y;
  const TVector3 yp = -sin(a) * x + cos(a) * y;
  const TVector3 zp = z.Unit();

  return {xp, yp, zp};
}

//____________________________________________________________________________
void EventLibraryInterface::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  LoadRecords();
}

//___________________________________________________________________________
void EventLibraryInterface::Configure(string config)
{
  Algorithm::Configure(config);
  LoadRecords();
}

//___________________________________________________________________________
void EventLibraryInterface::Cleanup()
{
  for(auto it: fRecords) delete it.second;
  fRecords.clear();
  delete fRecordFile;
  fRecordFile = 0;
}

//___________________________________________________________________________
void EventLibraryInterface::LoadRecords()
{
  Cleanup();

  std::string libPath;
  GetParam("EventLibraryPath", libPath);
  Expand(libPath);

  bool onDemand;
  GetParam("OnDemand", onDemand);

  PDGLibrary* pdglib = PDGLibrary::Instance();

  fRecordFile = new TFile(libPath.c_str());
  if(fRecordFile->IsZombie()) exit(1);

  TIter next(fRecordFile->GetListOfKeys());
  while(TObject* dir = next()){
    const std::string& tgtName = dir->GetName();
    const TParticlePDG* tgtPart = pdglib->DBase()->GetParticle(tgtName.c_str());
    if(!tgtPart){
      LOG("ELI", pWARN) << "Unknown nucleus " << tgtName
                        << " found in " << libPath
                        << " -- skipping";
      continue;
    }

    for(int sign: {+1, -1}){
      for(int pdg: {12, 14, 16}){
        for(bool iscc: {true, false}){
          // NCs should be the same for all flavours. Use nu_mu as a
          // convention internal to this code to index into the records map.
          if(!iscc && pdg != 14) continue;

          std::string nuName = pdglib->Find(sign*pdg)->GetName();
          if(!iscc) nuName = (sign > 0) ? "nu" : "nu_bar";

          const std::string treeName =
            TString::Format("%s/%s/%s/records",
                            tgtName.c_str(),
                            iscc ? "cc" : "nc",
                            nuName.c_str()).Data();

          const Key key(tgtPart->PdgCode(), sign*pdg, iscc);

          TTree* tr = (TTree*)fRecordFile->Get(treeName.c_str());

          if(!tr){
            LOG("ELI", pINFO) << treeName << " not found in "
                              << libPath << " -- skipping";
            continue;
          }

          if(onDemand)
            fRecords[key] = new OnDemandRecordList(tr, treeName);
          else
            fRecords[key] = new SimpleRecordList(tr, treeName);
        } // end for nucleus
      } // end for pdg
    } // end for sign
  } // end for iscc

  // Need to keep the record file open for OnDemand, but not Simple
  if(!onDemand){delete fRecordFile; fRecordFile = 0;}
}
