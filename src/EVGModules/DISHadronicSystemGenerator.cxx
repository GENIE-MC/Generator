//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMCParticle6.h>

#include "Conventions/Constants.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/DISHadronicSystemGenerator.h"
#include "Fragmentation/HadronizationModelI.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/FragmRecUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::print;

//___________________________________________________________________________
DISHadronicSystemGenerator::DISHadronicSystemGenerator() :
HadronicSystemGenerator("genie::DISHadronicSystemGenerator")
{

}
//___________________________________________________________________________
DISHadronicSystemGenerator::DISHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::DISHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
DISHadronicSystemGenerator::~DISHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void DISHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  //-- If the struck nucleon was within a nucleus, then add the final state
  //   nucleus at the EventRecord
  this->AddTargetNucleusRemnant(evrec);

  //-- Add an entry for the DIS Pre-Fragm. Hadronic State
  this->AddFinalHadronicSyst(evrec);

  //-- Add the fragmentation products
  this->AddFragmentationProducts(evrec);
}
//___________________________________________________________________________
void DISHadronicSystemGenerator::AddFragmentationProducts(
                                                    GHepRecord * evrec) const
{
// Calls a hadronizer and adds the fragmentation products at the GHEP

  //-- Compute the hadronic system invariant mass
  TLorentzVector p4Had = this->Hadronic4pLAB(evrec);
  double W = p4Had.M();

  Interaction * interaction = evrec->Summary();
  interaction->KinePtr()->SetW(W);

  //-- Run the hadronization model and get the fragmentation products:
  //   A collection of ROOT TMCParticles (equiv. to a LUJETS record)

  TClonesArray * plist = fHadronizationModel->Hadronize(interaction);
  if(!plist) {
     LOG("DISHadronicVtx", pWARN) 
                  << "Got an empty particle list. Hadronizer failed!";
     LOG("DISHadronicVtx", pWARN) 
                    << "Quitting the current event generation thread";

     evrec->EventFlags()->SetBitNumber(kNoAvailablePhaseSpace, true);

     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Not enough phase space for hadronizer");
     exception.SwitchOnFastForward();
     throw exception;

     return;
  }

  //-- Take the hadronic system weight to handle cases that the hadronizer
  //   was asked to produce weighted events
  double wght = fHadronizationModel->Weight();

  //-- Velocity for the [Hadronic CM] -> [LAB] active Lorentz transform
  TVector3 beta = this->HCM2LAB(evrec);

  //-- Translate the fragmentation products from TMCParticles to
  //   GHepParticles and copy them to the event record.

  int mom = evrec->FinalStateHadronicSystemPosition();
  assert(mom!=-1);

  TLorentzVector v4(0,0,0,0); // dummy position 4-vector
 
  TMCParticle * p = 0;
  TIter particle_iter(plist);

  bool is_nucleus = interaction->InitState().Tgt().IsNucleus();
  GHepStatus_t istfin = (is_nucleus) ?       
                 kIStHadronInTheNucleus : kIStStableFinalState;

  //-- Get a unit momentum along the momentum transfer direction \vec{q}
  //   at the [Hadronic CM] 
  TLorentzVector pq4 = this->MomentumTransferLAB(evrec);
  pq4.Boost(-beta);
  TVector3 unitvq = pq4.Vect().Unit();

  while( (p = (TMCParticle *) particle_iter.Next()) ) {

     int pdgc = p->GetKF();
     int ks   = p->GetKS();

     if(fFilterPreFragmEntries && ks!=1) continue;

     // The fragmentation products are generated in the final state
     // hadronic CM Frame (with the z>0 axis being the \vec{q} direction).
     // For each particle returned by the hadronizer:
     // - rotate its 3-momentum so that z := \vec{q} at hadronic CM frame
     // - boost it back to LAB frame

     TVector3 p3(p->GetPx(), p->GetPy(), p->GetPz());
     p3.RotateUz(unitvq); 

     TLorentzVector p4(p3, p->GetEnergy());
     p4.Boost(beta); 

     // copy final state particles to the event record
     GHepStatus_t ist = (ks==1) ? istfin : kIStDISPreFragmHadronicState;

     int im  = mom + 1 + p->GetParent();
     int ifc = (p->GetFirstChild() == -1) ? -1 : mom + 1 + p->GetFirstChild();
     int ilc = (p->GetLastChild()  == -1) ? -1 : mom + 1 + p->GetLastChild();

     evrec->AddParticle(pdgc, ist, im,-1, ifc, ilc, p4,v4);

  } // fragmentation-products-iterator

  //-- Handle the case that the hadronizer produced weighted events and
  //   take into account that the current event might be already weighted
  evrec->SetWeight (wght * evrec->Weight());

  plist->Delete();
  delete plist;
}
//___________________________________________________________________________
void DISHadronicSystemGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISHadronicSystemGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISHadronicSystemGenerator::LoadConfig(void)
{
// Load sub-algorithms and config data to reduce the number of registry
// lookups

  fHadronizationModel = 0;

  //-- Get the requested hadronization model
  fHadronizationModel = dynamic_cast<const HadronizationModelI *> (
        this->SubAlg("hadronization-alg-name", "hadronization-param-set"));

  assert(fHadronizationModel);

  //-- flag to determine whether we copy all fragmentation record entries
  //   into the GHEP record or just the ones marked with kf=1
  fFilterPreFragmEntries = 
                     fConfig->GetBoolDef("filter-pre-fragm-entries",false);
}
//____________________________________________________________________________

