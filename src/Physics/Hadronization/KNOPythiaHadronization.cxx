//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 10, 2011
   Fixed a bug reported by Torben Ferber affecting the KNO -> PYTHIA
   model transition (the order was reversed!)

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TH1D.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Physics/Hadronization/KNOPythiaHadronization.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
KNOPythiaHadronization::KNOPythiaHadronization() :
EventRecordVisitorI("genie::KNOPythiaHadronization")
{

}
//____________________________________________________________________________
KNOPythiaHadronization::KNOPythiaHadronization(string config) :
EventRecordVisitorI("genie::KNOPythiaHadronization", config)
{

}
//____________________________________________________________________________
KNOPythiaHadronization::~KNOPythiaHadronization()
{

}
//____________________________________________________________________________
// void KNOPythiaHadronization::Initialize(void) const
// {
//
// }
//____________________________________________________________________________
void KNOPythiaHadronization::ProcessEventRecord(GHepRecord * event) const
{
// Generate the hadronic system using either the KNO-based or PYTHIA/JETSET
// hadronization models according to the specified transition scheme
  Interaction * interaction = event->Summary();


  //-- Select hadronizer
  const EventRecordVisitorI * hadronizer =
      this->SelectHadronizer(interaction);

  //-- Run the selected hadronizer
  hadronizer->ProcessEventRecord(event);

  // //-- Update the weight
  // fWeight = hadronizer->Weight();
}
//____________________________________________________________________________
/*
PDGCodeList * KNOPythiaHadronization::SelectParticles(
                                        const Interaction * interaction) const
{
  //-- Select hadronizer
  const EventRecordVisitorI * hadronizer = this->SelectHadronizer(interaction);

  //-- Run the selected hadronizer
  PDGCodeList * pdgv = hadronizer->SelectParticles(interaction);

  return pdgv;
}
//____________________________________________________________________________
TH1D * KNOPythiaHadronization::MultiplicityProb(
 		       const Interaction * interaction, Option_t * opt) const
{
  //-- Select hadronizer
  const EventRecordVisitorI * hadronizer = this->SelectHadronizer(interaction);

  //-- Run the selected hadronizer
  TH1D * mprob = hadronizer->MultiplicityProb(interaction,opt);

  return mprob;
}
//____________________________________________________________________________
double KNOPythiaHadronization::Weight(void) const
{
  return fWeight;
}
*/
//____________________________________________________________________________
const EventRecordVisitorI * KNOPythiaHadronization::SelectHadronizer(
                                        const Interaction * interaction) const
{
  const EventRecordVisitorI * hadronizer = 0;
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
      if(R>f) hadronizer = fKNOHadronizer;
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

   // Load the requested hadronizers
  fKNOHadronizer =
     dynamic_cast<const EventRecordVisitorI *> (
                              this->SubAlg("KNO-Hadronizer"));
  fPythiaHadronizer =
     dynamic_cast<const EventRecordVisitorI *> (
                              this->SubAlg("PYTHIA-Hadronizer"));

  assert(fKNOHadronizer && fPythiaHadronizer);

  // Get transition method
  fMethod = 2 ;
  GetParam( "TransMethod", fMethod, false ) ;


  // Get transition scheme specific config
  if(fMethod==2) {

	GetParam( "KNO2PYTHIA-Wmin", fWminTrWindow ) ;

	GetParam( "KNO2PYTHIA-Wmax", fWmaxTrWindow ) ;

  }
}
//____________________________________________________________________________
