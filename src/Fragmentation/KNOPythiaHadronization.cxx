//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
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

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Fragmentation/KNOPythiaHadronization.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"

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

  double W = interaction->GetKinematics().W();
  LOG("HybridHad", pINFO) << "W = " << W << " GeV";

  if(W <= kNucleonMass+kPionMass) {
     LOG("HybridHad", pWARN) 
        << "Low invariant mass, W = " << W << " GeV! Returning a null list";
     return 0;
  }

  TClonesArray * particle_list = 0;

  switch(fMethod) {

  // ** KNO-only
  case(0) :
    particle_list = fKNOHadronizer->Hadronize(interaction);
    fWeight =  fKNOHadronizer->Weight();
    break;

  // ** PYTHIA/JETSET-only
  case(1) :
    particle_list = fPythiaHadronizer->Hadronize(interaction);
    fWeight =  fPythiaHadronizer->Weight();
    break;

  // ** KNO-only           @ W < Wmin
  // ** PYTHIA/JETSET-only @ W > Wmax
  // ** Smooth linear transition in [Wmin,Wmax]
  case(2) :
    particle_list = this->LinearTransitionWindowMethod(interaction);
    break;

  default :
    LOG("HybridHad", pFATAL) 
                    << "Unspecified transition method: " << fMethod;
    exit(1);
  }

  return particle_list;
}
//____________________________________________________________________________
TClonesArray * KNOPythiaHadronization::LinearTransitionWindowMethod(
                                        const Interaction * interaction) const
{
  TClonesArray * particle_list = 0;

  RandomGen * rnd = RandomGen::Instance();

  //----- Init event weight (to be set if producing weighted events)
  fWeight = 1.;

  double W = interaction->GetKinematics().W();

  if(W <= fWminTrWindow) {
    // KNO-only
    particle_list = fKNOHadronizer->Hadronize(interaction);
    fWeight =  fKNOHadronizer->Weight();
  } 

  else if (W > fWminTrWindow) {
    // PYTHIA/JETSET-only
    particle_list = fPythiaHadronizer->Hadronize(interaction);
    fWeight =  fPythiaHadronizer->Weight();

  } else {
    // Transition window
    double R = rnd->RndHadro().Rndm();
    double f = (W-fWminTrWindow)/(fWmaxTrWindow-fWminTrWindow);
    if(R<f) {
       particle_list = fKNOHadronizer->Hadronize(interaction);
       fWeight =  fKNOHadronizer->Weight();
    } else {
       particle_list = fPythiaHadronizer->Hadronize(interaction);
       fWeight =  fPythiaHadronizer->Weight();
    }
  }
  return particle_list;
}
//____________________________________________________________________________
double KNOPythiaHadronization::Weight(void) const
{
  return fWeight;
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

  // Get config sets for KNO and PYTHIA hadonizers
  string kno_config    = fConfig->GetString("kno-model-param-set");
  string pythia_config = fConfig->GetString("pythia-model-param-set");

  // Load the requested hadronizers
  AlgFactory * algf = AlgFactory::Instance();

  fKNOHadronizer = dynamic_cast<const HadronizationModelI *> (
                  algf->GetAlgorithm("genie::KNOHadronization", kno_config));
  fPythiaHadronizer = dynamic_cast<const HadronizationModelI *> (
            algf->GetAlgorithm("genie::PythiaHadronization", pythia_config));

  assert(fKNOHadronizer && fPythiaHadronizer);

  // Get transition method
  fMethod = fConfig->GetIntDef("transition-method", -1);

  // Get transition scheme specific config
  if(fMethod==2) {
    fWminTrWindow = fConfig->GetDouble("transition-window-Wmin");
    fWmaxTrWindow = fConfig->GetDouble("transition-window-Wmax");
  }
}
//____________________________________________________________________________
