//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 08, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TVector3.h>

#include "EVGCore/EVGThreadException.h"
#include "EVGModules/PauliBlocker.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepFlags.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Nuclear/FermiMomentumTablePool.h"
#include "Nuclear/FermiMomentumTable.h"

using namespace genie;

//___________________________________________________________________________
PauliBlocker::PauliBlocker() :
EventRecordVisitorI("genie::PauliBlocker")
{

}
//___________________________________________________________________________
PauliBlocker::PauliBlocker(string config) :
EventRecordVisitorI("genie::PauliBlocker",  config)
{

}
//___________________________________________________________________________
PauliBlocker::~PauliBlocker()
{

}
//___________________________________________________________________________
void PauliBlocker::ProcessEventRecord(GHepRecord * event_rec) const
{
  //-- Get the Interaction & InitialState objects
  Interaction * interaction = event_rec->Summary();
  const InitialState & init_state = interaction->InitState();

  //-- Pauli Blocking is only relevant for nucleon bound in a nucleus

  if( init_state.Tgt().IsNucleus() ) {

    int tgt_pdgc = init_state.Tgt().Pdg();
    int nuc_pdgc = interaction->RecoilNucleonPdg();

    GHepParticle * hit = event_rec->HitNucleon();
    TVector3 beta = hit->P4()->BoostVector();
  
    if(nuc_pdgc != 0) {
       // Find the recoil nucleon in the EventRecord

       GHepStatus_t ist = kIStStableFinalState;
       GHepParticle * nuc = event_rec->FindParticle(nuc_pdgc, ist, 0);

       if(nuc) {
         // get the Fermi momentum
         const double kf = fKFTable->FindClosestKF(tgt_pdgc, nuc_pdgc);

         LOG("PauliBlock", pINFO) << "KF = " << kf;

         //TLorentzVector * p4 = nuc->GetP4();
         //p4->Boost(-beta);
         //double p = p4->P(); // |p| for the recoil nucleon         
         //delete p4;

	 double p = nuc->P4()->P(); // |p| for the recoil nucleon
         LOG("PauliBlock", pINFO) << "Recoil nucleon |P| = " << p;

         if(p < kf) {
              LOG("PauliBlock", pINFO)
                   << "\n The generated event is Pauli-blocked: "
                          << " |p| = " << p << " < Fermi-Momentum = " << kf;

              const ProcessInfo & proc = interaction->ProcInfo();

              event_rec->EventFlags()->SetBitNumber(kPauliBlock, true);
              genie::exceptions::EVGThreadException exception;
              exception.SetReason("Pauli-blocked event");

              if(proc.IsQuasiElastic()) {
                 // nuclear suppression taken into account at the QEL cross
                 // section - should attempt to regenerate the event as QEL
                 exception.SwitchOnStepBack();
                 exception.SetReturnStep(0);
              } else {
                 // end this event generation thread and start again at the 
                 // interaction selection step
                 exception.SwitchOnFastForward();
              }
              throw exception;
         }
       }//nuc!=0
    }//nuc_pdgc!=0
  }//not a free nucleon
}
//___________________________________________________________________________
void PauliBlocker::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadKFTable();
}
//___________________________________________________________________________
void PauliBlocker::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadKFTable();
}
//___________________________________________________________________________
void PauliBlocker::LoadKFTable(void)
{
  fKFTable = 0;

  // get the Fermi momentum table
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  fKFTable = kftp->GetTable("Default");

  assert(fKFTable);
}
//___________________________________________________________________________

