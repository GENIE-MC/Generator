//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TVector3.h>

#include "Algorithm/AlgConfigPool.h"
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
void PauliBlocker::ProcessEventRecord(GHepRecord * evrec) const
{
  // Return if the neutrino was not scatterred off a nuclear target
  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("PauliBlock", pINFO)
           << "No nuclear target found - The Pauli Blocker exits";
    return;
  }

  // Handle only QEL for now...
  Interaction * interaction = evrec->Summary();
  const ProcessInfo & proc = interaction->ProcInfo();
  if(!proc.IsQuasiElastic()) {
    LOG("PauliBlock", pINFO) << "Not a QEL event - The Pauli Blocker exits";  
    return;
  }

  GHepParticle * hit = evrec->HitNucleon();
  assert(hit);
  GHepParticle * recoil = evrec->Particle(hit->FirstDaughter());
  assert(recoil);

  int tgt_pdgc = nucltgt -> Pdg();
  int nuc_pdgc = recoil  -> Pdg();

  // get the Fermi momentum
  const double kf = fKFTable->FindClosestKF(tgt_pdgc, nuc_pdgc);
  LOG("PauliBlock", pINFO) << "KF = " << kf;

  // get the recoil momentum
  double p = recoil->P4()->P(); // |p| for the recoil nucleon
  LOG("PauliBlock", pINFO) << "Recoil nucleon |P| = " << p;

  // check for pauli blocking
  bool is_blocked = (p < kf);

  // if it is blocked, set & thow an exception
  if(is_blocked) {
     LOG("PauliBlock", pNOTICE)
        << " *** The generated event is Pauli-blocked ("
        << "|p_{nucleon}| = " << p << " GeV < Fermi momentum = " << kf << " GeV) ***";

     evrec->EventFlags()->SetBitNumber(kPauliBlock, true);
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
        // - this is irrelevant for the time being as we only handle QEL-
        exception.SwitchOnFastForward();
     }
     throw exception;
  }
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
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fKFTableName = fConfig->GetStringDef ("FermiMomentumTable",
                                    gc->GetString("FermiMomentumTable"));
  fKFTable = 0;

  // get the Fermi momentum table
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  fKFTable = kftp->GetTable(fKFTableName);
  assert(fKFTable);
}
//___________________________________________________________________________

