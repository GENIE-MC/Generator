//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - August 17, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TH1D.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Fragmentation/KNOPythiaHadronization.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodeList.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
KNOPythiaHadronization::KNOPythiaHadronization() :
HadronizationModelI("genie::KNOPythiaHadronization")
{

}
//____________________________________________________________________________
KNOPythiaHadronization::KNOPythiaHadronization(string config) :
HadronizationModelI("genie::KNOPythiaHadronization", config)
{

}
//____________________________________________________________________________
KNOPythiaHadronization::~KNOPythiaHadronization()
{

}
//____________________________________________________________________________
void KNOPythiaHadronization::Initialize(void) const
{

}
//____________________________________________________________________________
TClonesArray * KNOPythiaHadronization::Hadronize(
                                        const Interaction * interaction) const
{
// Generate the hadronic system using either the KNO-based or PYTHIA/JETSET 
// hadronization models according to the specified transition scheme

  double W = interaction->Kine().W();
  LOG("HybridHad", pINFO) << "W = " << W << " GeV";

  if(W <= kNucleonMass+kPionMass) {
     LOG("HybridHad", pWARN) 
        << "Low invariant mass, W = " << W << " GeV! Returning a null list";
     return 0;
  }

  //-- Init event weight (to be set if producing weighted events)
  fWeight = 1.;

  //-- Select hadronizer
  const HadronizationModelI * hadronizer = this->SelectHadronizer(interaction);

  //-- Run the selected hadronizer
  TClonesArray * particle_list = hadronizer->Hadronize(interaction);

  //-- Update the weight
  fWeight = hadronizer->Weight();

  return particle_list;
}
//____________________________________________________________________________
PDGCodeList * KNOPythiaHadronization::SelectParticles(
                                        const Interaction * interaction) const
{
  //-- Select hadronizer
  const HadronizationModelI * hadronizer = this->SelectHadronizer(interaction);

  //-- Run the selected hadronizer
  PDGCodeList * pdgv = hadronizer->SelectParticles(interaction);

  return pdgv;
}
//____________________________________________________________________________
TH1D * KNOPythiaHadronization::MultiplicityProb(
 		       const Interaction * interaction, Option_t * opt) const
{
  //-- Select hadronizer
  const HadronizationModelI * hadronizer = this->SelectHadronizer(interaction);

  //-- Run the selected hadronizer
  TH1D * mprob = hadronizer->MultiplicityProb(interaction,opt);

  return mprob;
}
//____________________________________________________________________________
double KNOPythiaHadronization::Weight(void) const
{
  return fWeight;
}
//____________________________________________________________________________
const HadronizationModelI * KNOPythiaHadronization::SelectHadronizer(
                                        const Interaction * interaction) const
{
  const HadronizationModelI * hadronizer = 0;
  RandomGen * rnd = RandomGen::Instance();
  double W = 0;

  switch(fMethod) {

  // ** KNO-only
  case(0) :
    hadronizer = fKNOHadronizer;
    break;

  // ** PYTHIA/JETSET-only
  case(1) :
    hadronizer = fPythiaHadronizer;
    break;

  // ** KNO-only           @ W < Wmin
  // ** PYTHIA/JETSET-only @ W > Wmax
  // ** Smooth linear transition in [Wmin,Wmax]
  case(2) :
    W = interaction->Kine().W();
    if      (W <= fWminTrWindow) hadronizer = fKNOHadronizer;
    else if (W >  fWmaxTrWindow) hadronizer = fPythiaHadronizer;
    else {
      // Transition window
      double R = rnd->RndHadro().Rndm();
      double f = (W-fWminTrWindow)/(fWmaxTrWindow-fWminTrWindow);
      if(R<f) hadronizer = fKNOHadronizer;
      else    hadronizer = fPythiaHadronizer;
    }
    break;

  default :
    LOG("HybridHad", pFATAL) 
                    << "Unspecified transition method: " << fMethod;
    exit(1);
  }

  if(!hadronizer) {
    LOG("HybridHad", pFATAL) << "Null hadronizer!!";
    exit(1);
  }

  LOG("HybridHad", pINFO) << "Selected hadronizer: " << hadronizer->Id();
  return hadronizer;
}
//____________________________________________________________________________
void KNOPythiaHadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void KNOPythiaHadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void KNOPythiaHadronization::LoadConfig(void)
{
// Read configuration options or set defaults

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // Load the requested hadronizers
  fKNOHadronizer = 
     dynamic_cast<const HadronizationModelI *> (
                              this->SubAlg("KNO-Hadronizer"));
  fPythiaHadronizer = 
     dynamic_cast<const HadronizationModelI *> (
                              this->SubAlg("PYTHIA-Hadronizer"));

  assert(fKNOHadronizer && fPythiaHadronizer);

  // Get transition method
  fMethod = fConfig->GetIntDef("TransMethod", 2);

  // Get transition scheme specific config
  if(fMethod==2) {
    fWminTrWindow = fConfig->GetDoubleDef(
                "KNO2PYTHIA-Wmin", gc->GetDouble("KNO2PYTHIA-Wmin"));
    fWmaxTrWindow = fConfig->GetDoubleDef(
                "KNO2PYTHIA-Wmax", gc->GetDouble("KNO2PYTHIA-Wmax"));
  }
}
//____________________________________________________________________________
