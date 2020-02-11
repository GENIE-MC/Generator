//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GMode.h"
#include "Framework/Interaction/Target.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/NeutralHeavyLepton/NHLPrimaryVtxGenerator.h"
#include "Physics/NeutralHeavyLepton/NHLDecayUtils.h"
#include "Physics/NeutralHeavyLepton/NHLDecayMode.h"

using namespace genie;

//____________________________________________________________________________
NHLPrimaryVtxGenerator::NHLPrimaryVtxGenerator() :
EventRecordVisitorI("genie::NHLPrimaryVtxGenerator")
{

}
//____________________________________________________________________________
NHLPrimaryVtxGenerator::NHLPrimaryVtxGenerator(string config) :
EventRecordVisitorI("genie::NHLPrimaryVtxGenerator",config)
{

}
//____________________________________________________________________________
NHLPrimaryVtxGenerator::~NHLPrimaryVtxGenerator()
{

}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::ProcessEventRecord(GHepRecord * event) const
{
  Interaction * interaction = event->Summary();

  fCurrDecayMode = (NHLDecayMode_t) interaction->ExclTag().DecayMode();

  LOG("NHL", pNOTICE)
    << "Simulating NHL decay " << utils::nhl::AsString(fCurrDecayMode);

  this->AddInitialState(event);
  this->GenerateDecayProducts(event);
}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::AddInitialState(GHepRecord * event) const
{
  TLorentzVector v4(0,0,0,0);

  Interaction * interaction = event->Summary();
  double E = interaction->InitState().ProbeE(kRfLab);
  double M = PDGLibrary::Instance()->Find(kPdgNHL)->Mass();
  double p = TMath::Sqrt(E*E-M*M);

  TLorentzVector p4(0,0,p,E);

  event->AddParticle(kPdgNHL, kIStInitialState, 0,-1,-1,-1, p4, v4);
}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::GenerateDecayProducts(GHepRecord * event) const
{
  LOG("NHL", pINFO) << "Generating decay...";

  PDGCodeList pdgv = utils::nhl::DecayProductList(fCurrDecayMode);
  LOG("NHL", pINFO) << "Decay product IDs: " << pdgv;
  assert ( pdgv.size() >  1);

  LOG("NHL", pINFO) << "Performing a phase space decay...";

  // Get the decay product masses

  vector<int>::const_iterator pdg_iter;
  int i = 0;
  double * mass = new double[pdgv.size()];
  double   sum  = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
    int pdgc = *pdg_iter;
    double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
    mass[i++] = m;
    sum += m;
  }

  LOG("NHL", pINFO)
    << "Decaying N = " << pdgv.size() << " particles / total mass = " << sum;

  int nhl_id = 0;
  GHepParticle * nhl = event->Particle(nhl_id);
  assert(nhl);
  TLorentzVector * p4d = nhl->GetP4();
  TLorentzVector * v4d = nhl->GetX4();

  LOG("NHL", pINFO)
    << "Decaying system p4 = " << utils::print::P4AsString(p4d);

  // Set the decay
  bool permitted = fPhaseSpaceGenerator.SetDecay(*p4d, pdgv.size(), mass);
  if(!permitted) {
     LOG("NHL", pERROR)
       << " *** Phase space decay is not permitted \n"
       << " Total particle mass = " << sum << "\n"
       << " Decaying system p4 = " << utils::print::P4AsString(p4d);
     // clean-up
     delete [] mass;
     delete p4d;
     delete v4d;
     // throw exception
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Decay not permitted kinematically");
     exception.SwitchOnFastForward();
     throw exception;
  }

  // Get the maximum weight
  //double wmax = fPhaseSpaceGenerator.GetWtMax();
  double wmax = -1;
  for(int idec=0; idec<200; idec++) {
     double w = fPhaseSpaceGenerator.Generate();
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  wmax *= 2;

  LOG("NHL", pNOTICE)
     << "Max phase space gen. weight @ current hadronic system: " << wmax;

  // Generate an unweighted decay
  RandomGen * rnd = RandomGen::Instance();

  bool accept_decay=false;
  unsigned int itry=0;
  while(!accept_decay)
  {
     itry++;

     if(itry > controls::kMaxUnweightDecayIterations) {
       // report, clean-up and return
       LOG("NHL", pWARN)
           << "Couldn't generate an unweighted phase space decay after "
           << itry << " attempts";
       // clean up
       delete [] mass;
       delete p4d;
       delete v4d;
       // throw exception
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select decay after N attempts");
       exception.SwitchOnFastForward();
       throw exception;
     }
     double w  = fPhaseSpaceGenerator.Generate();
     if(w > wmax) {
        LOG("NHL", pWARN)
           << "Decay weight = " << w << " > max decay weight = " << wmax;
     }
     double gw = wmax * rnd->RndHadro().Rndm();
     accept_decay = (gw<=w);

     LOG("NHL", pINFO)
        << "Decay weight = " << w << " / R = " << gw
        << " - accepted: " << accept_decay;

  } //!accept_decay

  // Insert final state products into a TClonesArray of GHepParticle's
  TLorentzVector v4(*v4d);
  int idp = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
     int pdgc = *pdg_iter;
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(idp);
     GHepStatus_t ist = kIStStableFinalState;
     event->AddParticle(pdgc, ist, nhl_id,-1,-1,-1, *p4fin, v4);
     idp++;
  }

  // Clean-up
  delete [] mass;
  delete p4d;
  delete v4d;
}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void NHLPrimaryVtxGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void NHLPrimaryVtxGenerator::LoadConfig(void)
{

}
//___________________________________________________________________________
