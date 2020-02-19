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
#include "Physics/Hadronization/AGKY2019.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
AGKY2019::AGKY2019() :
EventRecordVisitorI("genie::AGKY2019")
{

}
//____________________________________________________________________________
AGKY2019::AGKY2019(string config) :
EventRecordVisitorI("genie::AGKY2019", config)
{

}
//____________________________________________________________________________
AGKY2019::~AGKY2019()
{

}
//____________________________________________________________________________
// void AGKY2019::Initialize(void) const
// {
//
// }
//____________________________________________________________________________
void AGKY2019::ProcessEventRecord(GHepRecord * event) const
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
PDGCodeList * AGKY2019::SelectParticles(
                                        const Interaction * interaction) const
{
  //-- Select hadronizer
  const EventRecordVisitorI * hadronizer = this->SelectHadronizer(interaction);

  //-- Run the selected hadronizer
  PDGCodeList * pdgv = hadronizer->SelectParticles(interaction);

  return pdgv;
}
//____________________________________________________________________________
TH1D * AGKY2019::MultiplicityProb(
 		       const Interaction * interaction, Option_t * opt) const
{
  //-- Select hadronizer
  const EventRecordVisitorI * hadronizer = this->SelectHadronizer(interaction);

  //-- Run the selected hadronizer
  TH1D * mprob = hadronizer->MultiplicityProb(interaction,opt);

  return mprob;
}
//____________________________________________________________________________
double AGKY2019::Weight(void) const
{
  return fWeight;
}
*/
//____________________________________________________________________________
const EventRecordVisitorI * AGKY2019::SelectHadronizer(
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
void AGKY2019::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AGKY2019::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AGKY2019::LoadConfig(void)
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
