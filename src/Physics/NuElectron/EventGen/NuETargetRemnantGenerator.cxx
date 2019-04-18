//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the NuE package from its previous location (EVGModules package)
 @ Feb 12, 2013 - CA (code from Rosen Matev)
   Handle the IMD annihilation channel.

*/
//____________________________________________________________________________

#include "Physics/NuElectron/EventGen/NuETargetRemnantGenerator.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
NuETargetRemnantGenerator::NuETargetRemnantGenerator() :
EventRecordVisitorI("genie::NuETargetRemnantGenerator")
{

}
//___________________________________________________________________________
NuETargetRemnantGenerator::NuETargetRemnantGenerator(string config) :
EventRecordVisitorI("genie::NuETargetRemnantGenerator", config)
{

}
//___________________________________________________________________________
NuETargetRemnantGenerator::~NuETargetRemnantGenerator()
{

}
//___________________________________________________________________________
void NuETargetRemnantGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  this -> AddElectronNeutrino     (evrec);
  this -> AddTargetNucleusRemnant (evrec);
}
//___________________________________________________________________________
void NuETargetRemnantGenerator::AddElectronNeutrino(GHepRecord * evrec) const
{
  //-- Get all initial & final state particles 4-momenta (in the LAB frame)

  //incoming v:
  GHepParticle * nu = evrec->Probe();

  //struck particle
  GHepParticle * el = evrec->HitElectron();

  //final state primary lepton:
  GHepParticle * l = evrec->FinalStatePrimaryLepton();

  assert(nu);
  assert(el);
  assert(l);

  //-- Force energy conservation
  //   Pv(Ev,pxv,pyv,pzv) + Pe(En,pxn,pyn,pzn) = Pfsl(El,pxl,pyl,pzl) + Px

  const TLorentzVector & p4v = *(nu->P4());
  const TLorentzVector & p4e = *(el->P4());
  const TLorentzVector & p4l = *(l->P4());
  const TLorentzVector & p4  = p4v + p4e - p4l; 

  //-- Vtx position
  const TLorentzVector & vtx = *(nu->X4());

  LOG("NuETargetRemnant", pINFO) << "Adding final state lepton from e- vtx";

  const ProcessInfo & proc_info = evrec->Summary()->ProcInfo();
  int mom  = evrec->HitElectronPosition();
  int pdgc = 0;
  if      (proc_info.IsNuElectronElastic()) pdgc = evrec->Summary()->InitState().ProbePdg();
  else if (proc_info.IsInverseMuDecay()) pdgc = kPdgNuE;
  else if (proc_info.IsIMDAnnihilation()) pdgc = kPdgAntiNuMu;
  assert(pdgc!=0);
  evrec->AddParticle(
           pdgc,kIStStableFinalState, mom,-1,-1,-1, p4, vtx);
}
//___________________________________________________________________________
void NuETargetRemnantGenerator::AddTargetNucleusRemnant(
                                                    GHepRecord * evrec) const
{
// add the remnant nuclear target at the GHEP record

  LOG("NuETargetRemnant", pDEBUG) << "Adding final state nucleus";

  //-- get A,Z for initial state nucleus
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  bool is_nucleus = init_state.Tgt().IsNucleus();
  if (!is_nucleus) {
    LOG("NuETargetRemnant", pDEBUG)
               << "Initial state not a nucleus - no remnant nucleus to add";
    return;
  }
  int A = init_state.Tgt().A();
  int Z = init_state.Tgt().Z();

  int    ipdgc = pdg::IonPdgCode(A, Z);
  double mass  = PDGLibrary::Instance()->Find(ipdgc)->Mass();

  //-- Add the nucleus to the event record
  LOG("NuETargetRemnant", pINFO)
          << "Adding nucleus [A = " << A << ", Z = " << Z
                                            << ", pdgc = " << ipdgc << "]";

  int mom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
             ipdgc,kIStStableFinalState, mom,-1,-1,-1, 0,0,0,mass, 0,0,0,0);
}
//___________________________________________________________________________
