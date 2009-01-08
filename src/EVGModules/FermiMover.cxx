//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 08, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "EVGModules/FermiMover.h"
#include "EVGCore/EVGThreadException.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepFlags.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearModelI.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
FermiMover::FermiMover() :
EventRecordVisitorI("genie::FermiMover")
{

}
//___________________________________________________________________________
FermiMover::FermiMover(string config) :
EventRecordVisitorI("genie::FermiMover", config)
{

}
//___________________________________________________________________________
FermiMover::~FermiMover()
{

}
//___________________________________________________________________________
void FermiMover::ProcessEventRecord(GHepRecord * event_rec) const
{
  Interaction *  interaction = event_rec   -> Summary();
  InitialState * init_state  = interaction -> InitStatePtr();
  Target *       tgt         = init_state  -> TgtPtr();

  // do nothing for non-nuclear targets
  if(!tgt->IsNucleus()) return;

  TLorentzVector * p4 = tgt->HitNucP4Ptr();

  // do nothing if the struct nucleon 4-momentum was set (eg as part of the
  // initial state selection)
  if(p4->Px()>0 || p4->Py()>0 || p4->Pz()>0) return;

  // generate a Fermi momentum & removal energy
  fNuclModel->GenerateNucleon(*tgt);
  TVector3 p3 = fNuclModel->Momentum3();
  double w    = fNuclModel->RemovalEnergy();

  LOG("FermiMover", pINFO) 
     << "Generated nucleon momentum: ("
     << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << ")";
  LOG("FermiMover", pINFO) 
     << "Generated nucleon removal energy: w = " << w;
  
  double pF2 = p3.Mag2(); // (fermi momentum)^2

  // access the hit nucleon and target nucleus at the GHEP record
  GHepParticle * nucleon = event_rec->HitNucleon();
  GHepParticle * nucleus = event_rec->TargetNucleus();
  assert(nucleon);
  assert(nucleus);

  nucleon->SetRemovalEnergy(w);

  // struck nucleon energy:
  // two possible prescriptions depending on whether you want to force
  // the sruck nucleon to be on the mass-shell or not...

  double EN=0;

  if(!fKeepNuclOnMassShell) {
     //-- compute A,Z for final state nucleus & get its PDG code 
     int nucleon_pdgc = nucleon->Pdg();
     bool is_p  = pdg::IsProton(nucleon_pdgc);
     int Z = (is_p) ? nucleus->Z()-1 : nucleus->Z();
     int A = nucleus->A() - 1;

     TParticlePDG * fnucleus = 0;
     int ipdgc = pdg::IonPdgCode(A, Z);
     fnucleus = PDGLibrary::Instance()->Find(ipdgc);
     if(!fnucleus) {
        LOG("FermiMover", pFATAL)
             << "No particle with [A = " << A << ", Z = " << Z
                            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
        exit(1);
     }
     //-- compute the energy of the struck (off the mass-shell) nucleus

     double Mf  = fnucleus -> Mass(); // remnant nucleus mass
     double Mi  = nucleus  -> Mass(); // initial nucleus mass

     EN = Mi - TMath::Sqrt(pF2 + Mf*Mf);

  } else {
     double MN  = nucleon->Mass();
     double MN2 = TMath::Power(MN,2);
     EN = TMath::Sqrt(MN2+pF2);
  }

  //-- update the struck nucleon 4p at the interaction summary and at
  //   the GHEP record
  p4->SetPx( p3.Px() );
  p4->SetPy( p3.Py() );
  p4->SetPz( p3.Pz() );
  p4->SetE ( EN      ); 

  nucleon->SetMomentum(*p4); // update GHEP value

  // sometimes, for interactions near threshold, Fermi momentum might bring
  // the neutrino energy in the nucleon rest frame below threshold (for the
  // selected interaction). In this case mark the event as unphysical and
  // abort the current thread.
  const KPhaseSpace & kps = interaction->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("FermiMover", pNOTICE)
                  << "Event below threshold after generating Fermi momentum";

     double Ethr = kps.Threshold();
     double Ev   = init_state->ProbeE(kRfHitNucRest);
     LOG("FermiMover", pNOTICE)
              << "Ev (@ nucleon rest frame) = " << Ev << ", Ethr = " << Ethr;

     event_rec->EventFlags()->SetBitNumber(kBelowThrNRF, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("E < Ethr after generating nucleon Fermi momentum");
     exception.SwitchOnFastForward();
     throw exception;
  }
}
//___________________________________________________________________________
void FermiMover::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FermiMover::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FermiMover::LoadConfig(void)
{
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed

//  AlgConfigPool * confp = AlgConfigPool::Instance();
//  const Registry * gc = confp->GlobalParameterList();

  fNuclModel = 0;

  RgKey nuclkey = "NuclearModel";
//  RgAlg nuclalg = fConfig->GetAlgDef(nuclkey, gc->GetAlg(nuclkey));
//  LOG("FermiMover", pINFO) << "Loading nuclear model: " << nuclalg;

  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);

  fKeepNuclOnMassShell = 
         fConfig->GetBoolDef("KeepHitNuclOnMassShell", false);
}
//____________________________________________________________________________

