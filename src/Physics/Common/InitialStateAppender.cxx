//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 14, 2009 - CA
   Include the struck e- for Glashow resonance reactions.
 @ Feb 12, 2013 - CA (code from Rosen Matev)
   Handle the IMD annihilation channel.

*/
//____________________________________________________________________________

#include <TLorentzVector.h>

#include "Physics/Common/InitialStateAppender.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
InitialStateAppender::InitialStateAppender() :
EventRecordVisitorI("genie::InitialStateAppender")
{

}
//___________________________________________________________________________
InitialStateAppender::InitialStateAppender(string config) :
EventRecordVisitorI("genie::InitialStateAppender", config)
{

}
//___________________________________________________________________________
InitialStateAppender::~InitialStateAppender()
{

}
//___________________________________________________________________________
void InitialStateAppender::ProcessEventRecord(GHepRecord * evrec) const
{
// Adds the initial state particles at the event record (the order is
// significant)

  LOG("ISApp", pINFO) << "Adding the initial state to the event record";

  //-- add the incoming neutrino to the event record
  this->AddNeutrino(evrec);

  //-- add the nuclear target at the event record (if any)
  this->AddNucleus(evrec);

  //-- add the struck nucleon to the event record (if any)
  //   It is added with status-code = 0 (init state) if the target was a
  //   free nucleon, or with a status-code = 11 (nucleon target) if the
  //   target was a nucleus.
  //   If the interaction was ve- elastic, inverse muon decay or Glashow 
  //   resonance then it will add the target e- instead.
  this->AddStruckParticle(evrec);
}
//___________________________________________________________________________
void InitialStateAppender::AddNeutrino(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfLab);
  const TLorentzVector v4(0.,0.,0.,0.);

  int pdgc = init_state.ProbePdg();

  LOG("ISApp", pINFO) << "Adding neutrino [pdgc = " << pdgc << "]";

  evrec->AddParticle(pdgc,kIStInitialState, -1,-1,-1,-1, *p4, v4);

  delete p4;
}
//___________________________________________________________________________
void InitialStateAppender::AddNucleus(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  bool is_nucleus = init_state.Tgt().IsNucleus();
  if(!is_nucleus) {
    LOG("ISApp", pINFO)
         << "Not an interaction with a nuclear target - no nucleus to add";
    return;
  }
  int    A    = init_state.Tgt().A();
  int    Z    = init_state.Tgt().Z();
  int    pdgc = pdg::IonPdgCode(A, Z);
  double M    = PDGLibrary::Instance()->Find(pdgc)->Mass();

  LOG("ISApp", pINFO)
          << "Adding nucleus [A = " << A << ", Z = " << Z
                                              << ", pdg = " << pdgc << "]";

  evrec->AddParticle(pdgc,kIStInitialState,-1,-1,-1,-1, 0,0,0,M, 0,0,0,0);
}
//___________________________________________________________________________
void InitialStateAppender::AddStruckParticle(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo & proc_info   = interaction->ProcInfo();

  bool hit_e = proc_info.IsInverseMuDecay()    ||
               proc_info.IsIMDAnnihilation()   ||
               proc_info.IsNuElectronElastic() ||
               proc_info.IsGlashowResonance();

  if(hit_e) {
    int    pdgc = kPdgElectron;
    double mass = PDGLibrary::Instance()->Find(pdgc)->Mass();
    const TLorentzVector p4(0,0,0, mass);
    const TLorentzVector v4(0.,0.,0.,0.);

    LOG("ISApp", pINFO) << "Adding struck electron";
    evrec->AddParticle(pdgc, kIStInitialState, 1, -1, -1, -1, p4, v4);
    return;
  }

  int pdgc = init_state.Tgt().HitNucPdg();

  if(pdgc != 0) {

    bool is_nucleus = init_state.Tgt().IsNucleus();

    GHepStatus_t ist   = (is_nucleus) ? kIStNucleonTarget : kIStInitialState;
    int          imom1 = (is_nucleus) ? 1 : -1;
    int          imom2 = -1;

    const TLorentzVector p4(init_state.Tgt().HitNucP4());
    const TLorentzVector v4(0.,0.,0.,0.);

    LOG("ISApp", pINFO)<< "Adding struck nucleon [pdgc = " << pdgc << "]";

    evrec->AddParticle(pdgc, ist, imom1, imom2, -1, -1, p4, v4);

  }//if struck nucleon was set
}
//___________________________________________________________________________
