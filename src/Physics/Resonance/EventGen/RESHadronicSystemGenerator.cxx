//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 07, 2009 - CA
   Removed call to AddTargetNucleusRemnant(). This simulation step is now
   performed further upstream in the processing chain.
 @ Mar 03, 2009 - CA
   Moved into the new RES package from its previous location (EVGModules)
 @ Jul 23, 2010 - CA
   Use ResonanceCharge() from base class. Function removed from utils::res.

*/
//____________________________________________________________________________

// #include <RVersion.h>
// #if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
// #include <TMCParticle.h>
// #else
// #include <TMCParticle6.h>
// #endif

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/SppChannel.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
// #include "Physics/Decay/DecayModelI.h"
#include "Physics/Resonance/EventGen/RESHadronicSystemGenerator.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
RESHadronicSystemGenerator::RESHadronicSystemGenerator() :
HadronicSystemGenerator("genie::RESHadronicSystemGenerator")
{

}
//___________________________________________________________________________
RESHadronicSystemGenerator::RESHadronicSystemGenerator(string config):
HadronicSystemGenerator("genie::RESHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
RESHadronicSystemGenerator::~RESHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void RESHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  // Get the right resonance PDG code so that the selected resonance
  // conserves charge
  int pdgc = GetResonancePdgCode(evrec);

  // Add the selected resonance
  this->AddResonance(evrec,pdgc);

  // Decay the resonance (and its decay products, if they include resonances)
  fResonanceDecayer->ProcessEventRecord(evrec);

  // Add the baryon resonance decay products at the event record
  //this->AddResonanceDecayProducts(evrec,pdgc);

  // Handle resonance decay channels to other resonances or short-living
  //   partices
  //LOG("RESHadronicVtx", pNOTICE)
  //   << "Decay any resonance in  the initial resonance decay products";
  //this->PreHadronTransportDecays(evrec);
}
//___________________________________________________________________________
int RESHadronicSystemGenerator::GetResonancePdgCode(GHepRecord * evrec) const
{
// In the RES thread the resonance is specifed when selecting interaction
// This method adds it to the GHEP record.

  Interaction * interaction = evrec->Summary();

  // Get resonance id
  const XclsTag & xcls = interaction->ExclTag();
  assert(xcls.KnownResonance());
  Resonance_t res = xcls.Resonance();

  // Get resonance charge
  int q_res = this->ResonanceCharge(evrec);

  // Find resonance PDG code from resonance charge and id
  int pdgc = utils::res::PdgCode(res, q_res);

  LOG("RESHadronicVtx", pNOTICE)
     << "Selected event has RES with PDGC = " << pdgc << ", Q = " << q_res;

  return pdgc;
}
//___________________________________________________________________________
void RESHadronicSystemGenerator::AddResonance(
                                         GHepRecord * evrec, int pdgc) const
{
  //-- Compute RES p4 = p4(neutrino) + p4(hit nucleon) - p4(primary lepton)
  TLorentzVector p4 = this->Hadronic4pLAB(evrec);

  //-- Add the resonance at the EventRecord
  GHepStatus_t ist = kIStPreDecayResonantState;
  int mom = evrec->HitNucleonPosition();

  //-- Get vtx position
  GHepParticle * neutrino  = evrec->Probe();
  const TLorentzVector & vtx = *(neutrino->X4());

  evrec->AddParticle(pdgc, ist, mom,-1,-1,-1, p4, vtx);
}
//___________________________________________________________________________
// void RESHadronicSystemGenerator::AddResonanceDecayProducts(
//    				         GHepRecord * evrec, int pdgc) const
// {
// // Decay the baryon resonance, take the decay products, boost them in the LAB
// // and add them in the GHEP record.
// // Unlike the SPP thread where the resonance decay products are determined
// // from the selected SPP channel, in the RES thread we can any of the the
// // resonance's kinematically available(the RES is not on the mass shell)decay
// // channels
//
//   // find the resonance position
//   int irpos = evrec->ParticlePosition(pdgc, kIStPreDecayResonantState, 0);
//   assert(irpos>0);
//
//   // access the GHEP entry
//   GHepParticle * resonance = evrec->Particle(irpos);
//   assert(resonance);
//
//   // resonance location
//   const TLorentzVector & x4 = *(resonance->X4());
//
//   // prepare the decayer inputs
//   DecayerInputs_t dinp;
//   dinp.PdgCode = pdgc;
//   dinp.P4      = resonance->P4();
//
//   // do the decay
//   TClonesArray * decay_products = fResonanceDecayer->Decay(dinp);
//   if(!decay_products) {
//      LOG("RESHadronicVtx", pWARN) << "Got an empty decay product list!";
//      LOG("RESHadronicVtx", pWARN)
//                       << "Quitting the current event generation thread";
//
//      evrec->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
//
//      genie::exceptions::EVGThreadException exception;
//      exception.SetReason("Not enough phase space for hadronizer");
//      exception.SwitchOnFastForward();
//      throw exception;
//
//      return;
//   }
//
//   // get the decay weight (if any)
//   double wght = fResonanceDecayer->Weight();
//
//   // update the event weight
//   evrec->SetWeight(wght * evrec->Weight());
//
//   // decide the istatus of decay products
//   GHepParticle * nuc = evrec->TargetNucleus();
//   GHepStatus_t dpist = (nuc) ? kIStHadronInTheNucleus : kIStStableFinalState;
//
//   // if the list is not empty, boost and copy the decay products in GHEP
//   if(decay_products) {
//
//      // first, mark the resonance as decayed
//      resonance->SetStatus(kIStDecayedState);
//
//      // loop over the daughter and add them to the event record
//      TMCParticle * dpmc = 0;
//      TObjArrayIter decay_iter(decay_products);
//
//      while( (dpmc = (TMCParticle *) decay_iter.Next()) ) {
//
//         int dppdg = dpmc->GetKF();
//         double px = dpmc->GetPx();
//         double py = dpmc->GetPy();
//         double pz = dpmc->GetPz();
//         double E  = dpmc->GetEnergy();
//         TLorentzVector p4(px,py,pz,E);
//
//        //-- Only add the decay products - the mother particle already exists
//        if(dpmc->GetKS()==1) {
//          evrec->AddParticle(dppdg,dpist,irpos,-1,-1,-1, p4, x4);
//        }
//      }
//
//      // done, release the original list
//      decay_products->Delete();
//      delete decay_products;
//   }// !=0
// }
//___________________________________________________________________________
void RESHadronicSystemGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void RESHadronicSystemGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void RESHadronicSystemGenerator::LoadConfig(void)
{
  fResonanceDecayer = 0;
  //fPreINukeDecayer  = 0;

  // Get the specified decayers
  fResonanceDecayer =
      dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("Decayer"));
  assert(fResonanceDecayer);
  // fPreINukeDecayer =
  //    dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("PreTransportDecayer"));
  // assert(fPreINukeDecayer);
}
//___________________________________________________________________________
