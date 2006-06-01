//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - November 23, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMCParticle6.h>

#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Decay/DecayModelI.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/RESHadronicSystemGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepFlags.h"
#include "Interaction/Interaction.h"
#include "Interaction/SppChannel.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"

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

  //-- Get the right resonance PDG code so that the selected resonance
  //   conserves charge
  int pdgc = GetResonancePdgCode(evrec);

  //-- Add the selected resonance
  this->AddResonance(evrec,pdgc);

  //-- If the struck nucleon was within a nucleus, then add the final state
  //   nucleus at the EventRecord
  this->AddTargetNucleusRemnant(evrec);

  //-- Add the baryon resonance decay products at the event record
  this->AddResonanceDecayProducts(evrec,pdgc);
}
//___________________________________________________________________________
int RESHadronicSystemGenerator::GetResonancePdgCode(GHepRecord * evrec) const
{
// In the RES thread the resonance is specifed when selecting interaction 
// This method adds it to the GHEP record.

  //-- Determine the RES pdg code (from the selected Resonance_t & charge)
  Interaction * interaction = evrec->GetInteraction();
  const XclsTag & xcls = interaction->GetExclusiveTag();
  assert(xcls.KnownResonance());
  Resonance_t res = xcls.Resonance();
  int charge = utils::res::ResonanceCharge(interaction);
  int pdgc   = utils::res::PdgCode(res,charge);

  LOG("RESHadronicVtx", pNOTICE)
      << "Selected event has RES with PDGC = " << pdgc << ", Q = " << charge;

  return pdgc;
}
//___________________________________________________________________________
void RESHadronicSystemGenerator::AddResonance(
                                         GHepRecord * evrec, int pdgc) const
{
  // compute RES p4 = p4(neutrino) + p4(hit nucleon) - p4(primary lepton)
  TLorentzVector p4 = this->Hadronic4pLAB(evrec);

  //-- Add the resonance at the EventRecord
  GHepStatus_t ist = kIStPreDecayResonantState;
  int mom = evrec->StruckNucleonPosition();

  evrec->AddParticle(
        pdgc, ist, mom,-1,-1,-1, p4.Px(),p4.Py(),p4.Pz(),p4.E(), 0,0,0,0);
}
//___________________________________________________________________________
void RESHadronicSystemGenerator::AddResonanceDecayProducts(
   				         GHepRecord * evrec, int pdgc) const
{
// Decay the baryon resonance, take the decay products, boost them in the LAB
// and add them in the GHEP record.
// Unlike the SPP thread where the resonance decay products are determined
// from the selected SPP channel, in the RES thread we can any of the the
// resonance's kinematically available(the RES is not on the mass shell)decay 
// channels

  // find the resonance position
  int irpos = evrec->ParticlePosition(pdgc, kIStPreDecayResonantState, 0);
  assert(irpos>0);

  // access the GHEP entry
  GHepParticle * resonance = evrec->Particle(irpos);
  assert(resonance);

  // prepare the decayer inputs
  DecayerInputs_t dinp;
  dinp.PdgCode = pdgc;
  dinp.P4      = resonance->P4();

  // do the decay
  TClonesArray * decay_products = fResonanceDecayer->Decay(dinp);
  if(!decay_products) {
     LOG("RESHadronicVtx", pWARN) << "Got an empty decay product list!";
     LOG("RESHadronicVtx", pWARN) 
                      << "Quitting the current event generation thread";

     evrec->EventFlags()->SetBitNumber(kNoAvailablePhaseSpace, true);

     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Not enough phase space for hadronizer");
     exception.SwitchOnFastForward();
     throw exception;

     return;
  }

  // get the decay weight (if any)
  double wght = fResonanceDecayer->Weight();

  // update the event weight
  evrec->SetWeight(wght * evrec->GetWeight());

  // decide the istatus of decay products
  GHepParticle * nuc = evrec->TargetNucleus();
  GHepStatus_t dpist = (nuc) ? kIStHadronInTheNucleus : kIStStableFinalState;

  // if the list is not empty, boost and copy the decay products in GHEP
  if(decay_products) {

     // first, mark the resonance as decayed
     resonance->SetStatus(kIStDecayedState);

     // loop over the daughter and add them to the event record
     TMCParticle * dpmc = 0;
     TObjArrayIter decay_iter(decay_products);

     while( (dpmc = (TMCParticle *) decay_iter.Next()) ) {

        int dppdg = dpmc->GetKF();
        double px = dpmc->GetPx();
        double py = dpmc->GetPy();
        double pz = dpmc->GetPz();
        double E  = dpmc->GetEnergy();

       //-- Only add the decay products - the mother particle already exists
       if(dpmc->GetKS()==1) {
         evrec->AddParticle(dppdg,dpist,irpos,-1,-1,-1, px,py,pz,E, 0,0,0,0);
       }
     }
         
     // done, release the original list
     decay_products->Delete();
     delete decay_products;
  }// !=0
}
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

  //-- Get the specified baryon resonance decayer
  fResonanceDecayer = dynamic_cast<const DecayModelI *>
                      (this->SubAlg("decayer-alg-name","decayer-param-set"));

  assert(fResonanceDecayer);
}
//___________________________________________________________________________

