//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org


*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
//#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
//#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Interaction/Interaction.h"
#include "Tools/VMC/EventLibraryInterface.h"
#include "Tools/VMC/RecordList.h"

#include <wordexp.h>

using namespace genie;
using namespace genie::vmc;

/// TODO there are almost certainly GENIE facilities for this
namespace
{
  const std::map<int, std::string> nucleus_to_label = {
  {1000010010, "H1"  },
  {1000060120, "C12" },
  {1000080160, "O16" },
  {1000170350, "Cl35"},
  {1000220480, "Ti48"},
  {1000260560, "Fe56"}};
}

//____________________________________________________________________________
EventLibraryInterface::EventLibraryInterface() :
EventRecordVisitorI("genie::vmc::EventLibraryInterface")
{

}
//____________________________________________________________________________
EventLibraryInterface::EventLibraryInterface(string config) :
EventRecordVisitorI("genie::EventLibraryInterface",config)
{

}
//____________________________________________________________________________
EventLibraryInterface::~EventLibraryInterface()
{
  for(auto it: fRecords) delete it.second;
}
//____________________________________________________________________________
void EventLibraryInterface::ProcessEventRecord(GHepRecord * event) const
{
// Get event summary constructed by GENIE
//
  Interaction * interaction = event->Summary();
  const InitialState & init_state = interaction->InitState();

  const Record* rec = GetRecord(init_state);
  if(!rec) return; // Reason has already been printed

  std::unique_ptr<TLorentzVector> probe_p4(init_state.GetProbeP4(kRfLab));

  event->AddParticle(init_state.ProbePdg(),
                     kIStInitialState,
                     -1,-1,-1,-1,
                     *probe_p4,
                     TLorentzVector(0, 0, 0, 0));

  const int tgt_A    = init_state.Tgt().A();
  const int tgt_Z    = init_state.Tgt().Z();
  const int tgt_pdgc = pdg::IonPdgCode(tgt_A, tgt_Z);
  // TODO why was that not just init_state.TgtPdg()?

  LOG("ELI", pINFO)
    << "Adding nucleus [A = " << tgt_A << ", Z = " << tgt_Z
    << ", pdg = " << tgt_pdgc << "]";

  event->AddParticle(tgt_pdgc,
                     kIStInitialState,
                     -1,-1,-1,-1,
                     0, 0, 0, PDGLibrary::Instance()->Find(tgt_pdgc)->Mass(),
                     0,0,0,0);

  const std::vector<TVector3> basis = Basis(probe_p4->Vect());

  for(const Particle& part: rec->parts){
    event->AddParticle(part.pdg,
                       kIStStableFinalState,
                       -1,-1,-1,-1,
                       TLorentzVector(part.px*basis[0] +
                                      part.py*basis[1] +
                                      part.pz*basis[2],
                                      part.E),
                       TLorentzVector(0, 0, 0, 0));
  }
}

//____________________________________________________________________________
const Record* EventLibraryInterface::GetRecord(const InitialState& init_state) const
{
  if(fRecords.empty()) LoadRecords();

  const std::unique_ptr<TLorentzVector> probe_p4(init_state.GetProbeP4(kRfLab));
  const double probe_E = probe_p4->E();
  const int probe_pdgc = init_state.ProbePdg();

  if(!init_state.Tgt().IsNucleus()){
    LOG("ELI", pINFO) << "Skippping non-nuclear target " << init_state;
    return 0;
  }

  const int tgt_A    = init_state.Tgt().A();
  const int tgt_Z    = init_state.Tgt().Z();
  const int tgt_pdgc = pdg::IonPdgCode(tgt_A, tgt_Z);

  const bool isCC = true; // TODO TODO this has already been decided, right? where do we find it?

  const Key key(tgt_pdgc, probe_pdgc, isCC);

  const auto rec_it = fRecords.find(key);

  if(rec_it == fRecords.end()){
    LOG("ELI", pINFO) << "Skippping " << key << " -- not found in library";
    return 0;
  }

  const Record* rec = rec_it->second->GetRecord(probe_E);

  if(!rec){
    LOG("ELI", pINFO) << "Skippping " << key << " at " << probe_E << " GeV -- not found in library";
    return 0;
  }

  return rec;
}

//____________________________________________________________________________
std::vector<TVector3> EventLibraryInterface::Basis(TVector3 z) const
{
  const TVector3 up(0, 1, 0);
  const TVector3 x = up.Cross(z).Unit(); // Perpendicular to neutrino and up
  const TVector3 y = x.Cross(z).Unit();  // Defines the third axis

  const double a = RandomGen::Instance()->RndEvg().Uniform(0, 2*M_PI);

  const TVector3 xp =  cos(a) * x + sin(a) * y;
  const TVector3 yp = -sin(a) * x + cos(a) * y;
  const TVector3 zp = z.Unit();

  return {xp, yp, zp};
}

//___________________________________________________________________________
// TODO is there a built-in tool for doing this?
void Expand(std::string& s)
{
  wordexp_t p;
  const int status = wordexp(s.c_str(), &p, WRDE_SHOWERR | WRDE_UNDEF);
  if(status != 0){
    LOG("ELI", pFATAL) << "String '" << s
                       << "' returned error " << status << " from wordexp().";
    exit(1);
  }

  if(p.we_wordc == 0){
    LOG("ELI", pFATAL) << "String '" << s
                       << "' didn't expand to anything.";
    exit(1);
  }

  if(p.we_wordc > 1){
    LOG("ELI", pFATAL) << "String '" << s
                       << "' expanded to " << p.we_wordc << " locations.";
    exit(1);
  }

  s = p.we_wordv[0];

  wordfree(&p);
}

//____________________________________________________________________________
void EventLibraryInterface::Configure(const Registry & config)
{
  Algorithm::Configure(config);
}
//___________________________________________________________________________
void EventLibraryInterface::Configure(string config)
{
  Algorithm::Configure(config);
}

//___________________________________________________________________________
void EventLibraryInterface::LoadRecords() const
{
  std::string dir;
  if(!GetParamDef( "LibraryPath", dir, std::string() )){
    LOG("ELI", pFATAL) << "Must specify 'LibraryPath'";
    exit(1);
  }

  Expand(dir);

  bool onDemand;
  GetParamDef( "OnDemand", onDemand, true );

  PDGLibrary* pdglib = PDGLibrary::Instance();

  for(bool iscc: {true, false}){
    for(int sign: {+1, -1}){
      for(int pdg: {12, 14}){
        // NCs should be the same for all flavours. Use numu by convention.
        if(!iscc && pdg == 12) continue;
        for(std::pair<int, std::string> it: nucleus_to_label){

            const TString fname =
              TString::Format("%s/%s%s_%s/%s/records.root",
                              dir.c_str(),
                              iscc ? "cc" : "nc",
                              (sign < 0) ? "bar" : "",
                              pdglib->Find(pdg)->GetName(),
                              it.second.c_str());

            const Key key(it.first, sign*pdg, iscc);

            if(onDemand)
              fRecords[key] = new OnDemandRecordList(fname.Data());
            else
              fRecords[key] = new SimpleRecordList(fname.Data());
        } // end for nucleus
      } // end for pdg
    } // end for sign
  } // end for iscc
}
