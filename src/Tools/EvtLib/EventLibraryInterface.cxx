//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org


*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Interaction/Interaction.h"
#include "Tools/EvtLib/EventLibraryInterface.h"
#include "Tools/EvtLib/EvtLibRecordList.h"
#include "Tools/EvtLib/Utils.h"
#include "Framework/Conventions/Constants.h"

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

  
  // at this point we have the event to add to our event record
  // the energy might be not the same
  // so first we correct the probe:
  // we maintain the same direction but we set the energy selected 
  // by the record.
  // We only change the the particle in the event record, 
  // Not the summary so that we keep the history
  // of the event construction
  
  // Neutrino is a parent to the lepton(s)
  GHepParticle * probe = event->Particle(0) ; 
  
  TLorentzVector probe_p4 ( * probe->P4() ) ;

  LOG("ELI", pINFO) << "Difference between requested neutrino E and used neutrino energy: " << probe_p4.E() - rec -> E ;

  probe_p4 *= rec -> E / probe_p4.E() ; 
  probe -> SetMomentum( probe_p4 ) ;

  // because of the selection procedure there might not be enough 
  // events in the library
  // so we add a roation around the z axis which is random for each event
  double alpha = RandomGen::Instance()->RndEvg().Uniform(0, 2*constants::kPi );

  
  // now we simply have to add the particles from the EvtLibRecord
  // But we need to rotate the particles along the direction of the beam
  // so we evaluate this direction once and for all
  TVector3 unit_nudir = probe->P4()->Vect().Unit();


  int firstLep =  -1 ; 
  const auto & parts = rec -> parts ;

  for ( unsigned int i = 0; i < parts.size() ; ++i ) {
    
    const auto & part = parts[i] ;
    
    int pdg = part.pdg;
    
    if(pdg::IsLepton(part.pdg)) {
      if ( firstLep == -1 ) { 
	firstLep = i ;
	interaction -> KinePtr() -> SetFSLeptonP4( part.px, part.py, part.pz, part.E ) ;
      }
    }
    
    // Fix up the outgoing lepton for NC events (due to lepton universality
    // it could be any value in the library)
    
    if( interaction->ProcInfo().IsWeakNC() ) {
      if ( int(i) == firstLep ) {
	// in the case of NC the process is always stored as numu (or numubar)
	// as the couplings are the same
	// here we restore the proper pdg on what we think it is the leading neutrino
	pdg = init_state.ProbePdg();
      }
    }

    TLorentzVector p4(part.px, part.py, part.pz, part.E) ;

    // The next is the rotation of the event according to the common alpha
    p4.RotateZ( alpha ) ;

    // Now we rotate the events so that it matches the direction of the incoming neutrino
    p4.RotateUz(unit_nudir); 

    event->AddParticle(pdg,
		       kIStStableFinalState,
		       0 , 1, -1, -1, // child of the neutrino and nucleus
		       p4,
		       TLorentzVector(0, 0, 0, 0));
    
  }
  
  // we now have the particle in the right direction inside the record we just need to 
  // decorate the summary kinematics
  
  if ( firstLep < 0 ) {
    firstLep = 2 ;
  }
  else { 
    firstLep += 2 ; 
  }
  // these are events from outside, there are no assumptions on the lepton
  // If we are simulating BSM physics there might not even be a lepton in the final state
  // So if there is a lepton we evaluate Q2 based on that lepton
  // If not we evaluate Q2 using the first particle we have available 
  
  FillKinematics( *event, *interaction->KinePtr(), firstLep );
  
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
    if(probe_pdgc > 0) probe_pdgc = kPdgNuMu; 
    else probe_pdgc = kPdgAntiNuMu;
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

    for(int pdg: {kPdgNuE,   kPdgAntiNuE,
                  kPdgNuMu,  kPdgAntiNuMu,
                  kPdgNuTau, kPdgAntiNuTau}){

      for(bool iscc: {true, false}){
        // NCs should be the same for all flavours. Use nu_mu as a convention
        // internal to this code to index into the records map.
        if(!iscc && abs(pdg) != kPdgNuMu) continue;

        std::string nuName = pdglib->Find(pdg)->GetName();
        if(!iscc) nuName = pdg::IsAntiNeutrino(pdg) ? "nu_bar" : "nu";

        const std::string treeName =
          TString::Format("%s/%s/%s/records",
                          tgtName.c_str(),
                          iscc ? "cc" : "nc",
                          nuName.c_str()).Data();

        const Key key(tgtPart->PdgCode(), pdg, iscc);

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
      } // end for iscc
    } // end for pdg
  } // end for dir

  // Need to keep the record file open for OnDemand, but not Simple
  if(!onDemand){delete fRecordFile; fRecordFile = 0;}
}

//___________________________________________________________________________
void EventLibraryInterface::FillKinematics( const GHepRecord & event,
					    Kinematics& kine, 
					    int primary_lep_id ) const {
  

  const TLorentzVector & p_probe = * event.Particle( 0 ) -> P4() ; 

  const TLorentzVector & p_lep = * event.Particle( primary_lep_id ) -> P4() ; 
  
  const TLorentzVector q = p_probe - p_lep;

  // Initial hadronic state, semi-arbitrary
  const TLorentzVector & p_tgt = * event.Particle( 1 ) -> P4() ;
  
  kine.SetQ2(-q.Mag2(), true);
  kine.SetW((q + p_tgt).Mag(), true);
  kine.Setx(-q.Mag2() / (2*p_tgt.Dot(q)), true);
  kine.Sety( p_tgt.Dot( q ) / p_tgt.Dot( p_probe ), true ) ;

}
