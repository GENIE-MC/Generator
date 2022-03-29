//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
 John Plows <komninos-john.plows \at physics.ox.ac.uk>
 University of Oxford
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
using namespace genie::NHL;

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

  fCurrInitStatePdg = interaction->InitState().ProbePdg();
  fCurrDecayMode = (NHLDecayMode_t) interaction->ExclTag().DecayMode();

  // interaction->ExclTag().SetHNL(); // need to modify Interaction/XclsTag

  LOG("NHL", pNOTICE)
    << "Simulating NHL decay " << utils::nhl::AsString(fCurrDecayMode)
    << " for an initial state with PDG code " << fCurrInitStatePdg;

  this->AddInitialState(event);
  //this->GenerateDecayPosition(event);
  this->GenerateDecayProducts(event);
}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::AddInitialState(GHepRecord * event) const
{
  TLorentzVector v4(0,0,0,0); // RETHERE make a way to sample the decay vertex position

  if( fUe42 == -1.0 && fUm42 == -1.0 && fUt42 == -1.0 ){
    LOG( "NHL", pINFO )
      << "Setting couplings to (1,1,0). This will change.";
    
    // SetNHLCouplings( event->Vertex()->X(), event->Vertex()->Y(), event->Vertex()->Z() );
    fUe42 = 1.0;
    fUm42 = 1.0;
    fUt42 = 0.0;
  }

  Interaction * interaction = event->Summary();
  double E = interaction->InitState().ProbeE(kRfLab);
  double M = PDGLibrary::Instance()->Find(kPdgNHL)->Mass();
  double p = TMath::Sqrt(E*E-M*M);

  // set some initial deviation from beam axis due to collimation effect
  // RETHERE make this configurable
  double thetaDev = 0.06; // deg
  thetaDev *= genie::constants::kPi / 180.0; // rad
  RandomGen * Rng = RandomGen::Instance();
  double theta = Rng->RndGen().Gaus(0.0, thetaDev);
  if( theta < 0.0 ) theta *= -1.0;
  double phi = Rng->RndGen().Uniform(0.0, 2.0 * genie::constants::kPi);

  double px = p * std::sin(theta) * std::cos(phi);
  double py = p * std::sin(theta) * std::sin(phi);
  double pz = p * std::cos(theta);
  
  TLorentzVector p4(px,py,pz,E);

  int hpdg = interaction->InitState().ProbePdg();
  event->AddParticle(hpdg, kIStInitialState, 0,-1,-1,-1, p4, v4);
}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::GenerateDecayProducts(GHepRecord * event) const
{
  LOG("NHL", pINFO) << "Generating decay...";

  // do we have nubar?
  PDGCodeList pdgv0 = utils::nhl::DecayProductList(fCurrDecayMode);
  int typeMod = ( fCurrInitStatePdg >= 0 ) ? 1 : -1; 
  PDGCodeList pdgv(true);
  for( std::vector<int>::iterator it = pdgv0.begin(); it != pdgv0.end(); ++it ){
    int pdgc = *it; 
    int newpdgc = ( pdgc == genie::kPdgPi0 ) ? pdgc : typeMod * pdgc; // pi-0 is its own antiparticle
    LOG("NHL", pDEBUG) << "Adding " << pdgc << " --> " << newpdgc;
    pdgv.push_back( newpdgc );
  }

  LOG("NHL", pINFO) << "Decay product IDs: " << pdgv;
  assert ( pdgv.size() > 1);

  // RETHERE may not want a phase space decay!
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
void NHLPrimaryVtxGenerator::SetNHLCouplings( double Ue42, double Um42, double Ut42 ) const
{
  fUe42 = Ue42;
  fUm42 = Um42;
  fUt42 = Ut42;
}
