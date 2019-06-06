//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Jeremy Hewes, Georgia Karagiorgi
         University of Manchester

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NNBarOscillation/NNBarOscPrimaryVtxGenerator.h"
#include "Physics/NNBarOscillation/NNBarOscUtils.h"
#include "Physics/NNBarOscillation/NNBarOscMode.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/NuclearState/NuclearModelI.h"

using namespace genie;

//____________________________________________________________________________
NNBarOscPrimaryVtxGenerator::NNBarOscPrimaryVtxGenerator() :
EventRecordVisitorI("genie::NNBarOscPrimaryVtxGenerator")
{

}
//____________________________________________________________________________
NNBarOscPrimaryVtxGenerator::NNBarOscPrimaryVtxGenerator(
  string config) :
EventRecordVisitorI("genie::NNBarOscPrimaryVtxGenerator",config)
{

}
//____________________________________________________________________________
NNBarOscPrimaryVtxGenerator::~NNBarOscPrimaryVtxGenerator() 
{

}
//____________________________________________________________________________
void NNBarOscPrimaryVtxGenerator::ProcessEventRecord(GHepRecord * event) const
{
  // spit out some output
  Interaction * interaction = event->Summary();
  fCurrInitStatePdg = interaction->InitState().Tgt().Pdg();
  fCurrDecayMode = (NNBarOscMode_t) interaction->ExclTag().DecayMode();

  // spit out that info -j
  LOG("NNBarOsc", pNOTICE)
    << "Simulating decay " << genie::utils::nnbar_osc::AsString(fCurrDecayMode)
    << " for an initial state with code: " << fCurrInitStatePdg;

  // check if nucleon is bound
  fNucleonIsBound = (pdg::IonPdgCodeToA(fCurrInitStatePdg) > 1);
  // can take this out for nnbar, nucleon is always bound!

  // now call these four functions!
  this->AddInitialState(event);
  this->GenerateOscillatingNeutronPosition(event);
  this->GenerateFermiMomentum(event);
  this->GenerateDecayProducts(event);
}
//____________________________________________________________________________
void NNBarOscPrimaryVtxGenerator::AddInitialState(GHepRecord * event) const
{

// add initial state to event record

// event record looks like this:
//    id: 0, nucleus (kIStInitialState)
//    |     
//    |---> { id: 1, neutron         (kIStDecayedState)
//          { id: 2, nucleon         (kIStDecayedState)
//          { id: 3, remnant nucleus (kIStStableFinalState)
//

  TLorentzVector v4(0,0,0,0);
  
  GHepStatus_t stis = kIStInitialState;
  GHepStatus_t stdc = kIStDecayedState;
  GHepStatus_t stfs = kIStStableFinalState;

  int ipdg = fCurrInitStatePdg;
  
  // add initial nucleus
  double Mi  = PDGLibrary::Instance()->Find(ipdg)->Mass();
  TLorentzVector p4i(0,0,0,Mi);
  event->AddParticle(ipdg,stis,-1,-1,-1,-1, p4i, v4);

  // add oscillating neutron
  int neutpdg = kPdgNeutron;
  double mneut = PDGLibrary::Instance()->Find(neutpdg)->Mass();
  TLorentzVector p4neut(0,0,0,mneut);
  event->AddParticle(neutpdg,stdc,0,-1,-1,-1, p4neut, v4);

  // add annihilation nucleon
  int dpdg = genie::utils::nnbar_osc::AnnihilatingNucleonPdgCode(fCurrDecayMode);
  double mn = PDGLibrary::Instance()->Find(dpdg)->Mass();
  TLorentzVector p4n(0,0,0,mn);
  event->AddParticle(dpdg,stdc, 0,-1,-1,-1, p4n, v4);

  // add nuclear remnant
  int A = pdg::IonPdgCodeToA(ipdg);
  int Z = pdg::IonPdgCodeToZ(ipdg);
  A--; A--;
  if(dpdg == kPdgProton) { Z--; }
  int rpdg = pdg::IonPdgCode(A, Z);
  double Mf  = PDGLibrary::Instance()->Find(rpdg)->Mass();
  TLorentzVector p4f(0,0,0,Mf);
  event->AddParticle(rpdg,stfs,0,-1,-1,-1, p4f, v4);
}
//____________________________________________________________________________
void NNBarOscPrimaryVtxGenerator::GenerateOscillatingNeutronPosition(
  GHepRecord * event) const
{
  GHepParticle * initial_nucleus = event->Particle(0);
  int A = initial_nucleus->A();
  if(A <= 2) {
    return;
  }

  RandomGen * rnd = RandomGen::Instance();

  double R0 = 1.3;
  double dA = (double)A;
  double R = R0 * TMath::Power(dA, 1./3.);
            
  LOG("NNBarOsc", pINFO)
      << "Generating vertex according to a realistic nuclear density profile";

  // get inputs to the rejection method
  double ymax = -1;
  double rmax = 3*R;
  double dr   = R/40.;
  for(double r = 0; r < rmax; r+=dr) {
      ymax = TMath::Max(ymax, r*r * utils::nuclear::Density(r,A));
  }
  ymax *= 1.2;
  
  // select a vertex using the rejection method 
  TLorentzVector vtx(0,0,0,0);
  unsigned int iter = 0;
  while(1) {
    iter++;

    // throw an exception if it hasn't find a solution after many attempts
    if(iter > controls::kRjMaxIterations) {
       LOG("NNBarOsc", pWARN)
           << "Couldn't generate a vertex position after " << iter << " iterations";
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't generate vertex");
       exception.SwitchOnFastForward();
       throw exception;
    }
           
    double r = rmax * rnd->RndFsi().Rndm();
    double t = ymax * rnd->RndFsi().Rndm();
    double y = r*r * utils::nuclear::Density(r,A);
    if(y > ymax) {   
       LOG("NNBarOsc", pERROR)
          << "y = " << y << " > ymax = " << ymax << " for r = " << r << ", A = " << A;
    }
    bool accept = (t < y);
    if(accept) {
       double phi      = 2*constants::kPi * rnd->RndFsi().Rndm();
       double cosphi   = TMath::Cos(phi);
       double sinphi   = TMath::Sin(phi);
       double costheta = -1 + 2 * rnd->RndFsi().Rndm();
       double sintheta = TMath::Sqrt(1-costheta*costheta);
       vtx.SetX(r*sintheta*cosphi);
       vtx.SetY(r*sintheta*sinphi);
       vtx.SetZ(r*costheta);
       vtx.SetT(0.);
       break;
    }
  } // while 1

  // giving position to oscillating neutron
  GHepParticle * oscillating_neutron = event->Particle(1);
  assert(oscillating_neutron);
  oscillating_neutron->SetPosition(vtx);

  // give same position to annihilation nucleon
  GHepParticle * annihilation_nucleon = event->Particle(2);
  assert(annihilation_nucleon);
  annihilation_nucleon->SetPosition(vtx);
}
//____________________________________________________________________________
void NNBarOscPrimaryVtxGenerator::GenerateFermiMomentum(
  GHepRecord * event) const
{
  GHepParticle * initial_nucleus = event->Particle(0);
  int A = initial_nucleus->A();
  if(A <= 2) {
    return;
  }

  GHepParticle * oscillating_neutron = event->Particle(1);
  GHepParticle * annihilation_nucleon = event->Particle(2);
  GHepParticle * remnant_nucleus = event->Particle(3);
  assert(oscillating_neutron);
  assert(annihilation_nucleon);
  assert(remnant_nucleus);

  // generate a Fermi momentum & removal energy
  Target tgt(initial_nucleus->Pdg());

  // start with oscillating neutron
  tgt.SetHitNucPdg(kPdgNeutron);
  // generate nuclear model & fermi momentum
  fNuclModel->GenerateNucleon(tgt);
  TVector3 p3 = fNuclModel->Momentum3();
  double w = fNuclModel->RemovalEnergy();

  // use mass & momentum to calculate energy
  double mass = oscillating_neutron->Mass();
  double energy = sqrt(pow(mass,2) + p3.Mag2()) - w;
  // give new energy & momentum to particle
  TLorentzVector p4(p3, energy);
  oscillating_neutron->SetMomentum(p4);

  LOG("NNBarOsc", pINFO)
     << "Generated neutron momentum: ("
     << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << "), "
     << "|p| = " << p3.Mag();

  // then rinse repeat for the annihilation nucleon
  tgt.SetHitNucPdg(annihilation_nucleon->Pdg());
  // use nuclear model to generate fermi momentum
  fNuclModel->GenerateNucleon(tgt);
  p3 = fNuclModel->Momentum3();
  w = fNuclModel->RemovalEnergy();
  // use mass & momentum to figure out energy
  mass = annihilation_nucleon->Mass();
  energy = sqrt(pow(mass,2) + p3.Mag2()) - w;
  // give new energy & momentum to particle
  p4 = TLorentzVector(p3, energy);
  annihilation_nucleon->SetMomentum(p4);

  LOG("NNBarOsc", pINFO) 
     << "Generated nucleon momentum: ("
     << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << "), "
     << "|p| = " << p3.Mag();

  // now figure out momentum for the nuclear remnant
  p3 = -1 * (oscillating_neutron->GetP4()->Vect() + annihilation_nucleon->GetP4()->Vect());
  // figure out energy from mass & momentum
  mass = remnant_nucleus->Mass();
  energy = sqrt(pow(mass,2) + p3.Mag2());
  // give new energy & momentum to remnant
  p4 = TLorentzVector(p3, energy);
  remnant_nucleus->SetMomentum(p4);
}
//____________________________________________________________________________
void NNBarOscPrimaryVtxGenerator::GenerateDecayProducts(
  GHepRecord * event) const
{
  LOG("NNBarOsc", pINFO) << "Generating decay...";

  PDGCodeList pdgv = genie::utils::nnbar_osc::DecayProductList(fCurrDecayMode);
  LOG("NNBarOsc", pINFO) << "Decay product IDs: " << pdgv;
  assert ( pdgv.size() >  1);

  LOG("NNBarOsc", pINFO) << "Performing a phase space decay...";

  // Get the decay product masses

  vector<int>::const_iterator pdg_iter;
  int idx = 0;
  double * mass = new double[pdgv.size()];
  double   sum  = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
    int pdgc = *pdg_iter;
    double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
    mass[idx++] = m;
    sum += m;
  }

  LOG("NNBarOsc", pINFO)  
    << "Decaying N = " << pdgv.size() << " particles / total mass = " << sum;
  int initial_nucleus_id      = 0;
  int oscillating_neutron_id  = 1;
  int annihilation_nucleon_id = 2;

  // get our annihilating nucleons
  GHepParticle * initial_nucleus      = event->Particle(initial_nucleus_id);
  assert(initial_nucleus);
  GHepParticle * oscillating_neutron  = event->Particle(oscillating_neutron_id);
  assert(oscillating_neutron);
  GHepParticle * annihilation_nucleon = event->Particle(annihilation_nucleon_id);
  assert(annihilation_nucleon);

  Target tgt(initial_nucleus->Pdg());
  tgt.SetHitNucPdg(kPdgNeutron);

  // get their momentum 4-vectors and boost into rest frame
  TLorentzVector * p4_1 = oscillating_neutron->GetP4();
  TLorentzVector * p4_2 = annihilation_nucleon->GetP4();
  TLorentzVector * p4d = new TLorentzVector(*p4_1 + *p4_2);
  TVector3 boost = p4d->BoostVector();
  p4d->Boost(-boost);

  // get decay position
  TLorentzVector * v4d = annihilation_nucleon->GetX4();

  delete p4_1;
  delete p4_2;

  LOG("NNBarOsc", pINFO) 
    << "Decaying system p4 = " << utils::print::P4AsString(p4d);

  // Set the decay
  bool permitted = fPhaseSpaceGenerator.SetDecay(*p4d, pdgv.size(), mass);

  // If the decay is not energetically allowed, select a new final state
  while(!permitted) {

    LOG("NNBarOsc", pINFO)
      << "Not enough energy to generate decay products! Selecting a new final state.";
    
    int mode = 0;
    
    int initial_nucleus_pdg = initial_nucleus->Pdg();
        
    std::string nucleus_pdg = std::to_string(static_cast<long long>(initial_nucleus_pdg));
    if (nucleus_pdg.size() != 10) {
      LOG("NNBarOsc", pERROR)
        << "Expecting the nuclear PDG code to be a 10-digit integer, but it is " << nucleus_pdg << ", which is "
        << nucleus_pdg.size() << " digits long. Drop me an email at jezhewes@gmail.com ; exiting...";
      exit(1);
    }

    int n_nucleons = std::stoi(nucleus_pdg.substr(6,3)) - 1;
    int n_protons  = std::stoi(nucleus_pdg.substr(3,3));
    double proton_frac  = ((double)n_protons) / ((double) n_nucleons);
    double neutron_frac = 1 - proton_frac;

    // set branching ratios, taken from bubble chamber data
    const int n_modes = 16;
    double br [n_modes] = { 0.010, 0.080, 0.100, 0.220,
                            0.360, 0.160, 0.070, 0.020,
                            0.015, 0.065, 0.110, 0.280,
                            0.070, 0.240, 0.100, 0.100 };
   
    for (int i = 0; i < n_modes; i++) {
      // for first 7 branching ratios, multiply by relative proton density
      if (i < 7)
        br[i] *= proton_frac;
      // for next 9, multiply by relative neutron density
      else
        br[i] *= neutron_frac;
    }

    // randomly generate a number between 1 and 0
    RandomGen * rnd = RandomGen::Instance();
    rnd->SetSeed(0);
    double p = rnd->RndNum().Rndm();
    
    // loop through all modes, figure out which one our random number corresponds to
    double threshold = 0;
    for (int j = 0; j < n_modes; j++) {
      threshold += br[j];
      if (p < threshold) {
        // once we've found our mode, stop looping
        mode = j + 1;
        break;
      }
    }
    
    // create new event record with new final state
    // TODO - I don't think Jeremy meant to make a _new_ record here; it
    // shadows the one passed in...
    // EventRecord * event = new EventRecord;
    Interaction * interaction = Interaction::NOsc((int)fCurrInitStatePdg, mode);
    event->AttachSummary(interaction);
    
    fCurrDecayMode = (NNBarOscMode_t) interaction->ExclTag().DecayMode(); 
    
    pdgv = genie::utils::nnbar_osc::DecayProductList(fCurrDecayMode);
    LOG("NNBarOsc", pINFO) << "Decay product IDs: " << pdgv;
    assert ( pdgv.size() > 1);
    
    // get the decay particles again
    LOG("NNBarOsc", pINFO) << "Performing a phase space decay...";
    idx = 0;
    delete [] mass;
    mass = new double[pdgv.size()];
    sum = 0;
    for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
      int pdgc = *pdg_iter;
      double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
      mass[idx++] = m;
      sum += m;
    }
    
    LOG("NNBarOsc", pINFO)
      << "Decaying N = " << pdgv.size() << " particles / total mass = " << sum;
    LOG("NNBarOsc", pINFO)
      << "Decaying system p4 = " << utils::print::P4AsString(p4d);
    
    permitted = fPhaseSpaceGenerator.SetDecay(*p4d, pdgv.size(), mass);
  }
  
  if(!permitted) {
     LOG("NNBarOsc", pERROR) 
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
  for(int i=0; i<200; i++) {
     double w = fPhaseSpaceGenerator.Generate();   
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  wmax *= 2;

  LOG("NNBarOsc", pNOTICE) 
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
       LOG("NNBarOsc", pWARN) 
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
        LOG("NNBarOsc", pWARN) 
           << "Decay weight = " << w << " > max decay weight = " << wmax;
     }
     double gw = wmax * rnd->RndHadro().Rndm();
     accept_decay = (gw<=w);

     LOG("NNBarOsc", pINFO) 
        << "Decay weight = " << w << " / R = " << gw 
        << " - accepted: " << accept_decay;

  } //!accept_decay

  // Insert final state products into a TClonesArray of TMCParticles
  TLorentzVector v4(*v4d); 
  int idp = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
     int pdgc = *pdg_iter;
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(idp);
     GHepStatus_t ist = 
        utils::nnbar_osc::DecayProductStatus(fNucleonIsBound, pdgc);
     p4fin->Boost(boost);
     event->AddParticle(pdgc, ist, oscillating_neutron_id,-1,-1,-1, *p4fin, v4);
     idp++;
  }

  // Clean-up
  delete [] mass;
  delete p4d;
  delete v4d;
}
//___________________________________________________________________________
void NNBarOscPrimaryVtxGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);   
  this->LoadConfig();
}
//___________________________________________________________________________
void NNBarOscPrimaryVtxGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void NNBarOscPrimaryVtxGenerator::LoadConfig(void)
{
//  AlgConfigPool * confp = AlgConfigPool::Instance();
//  const Registry * gc = confp->GlobalParameterList();
    
  fNuclModel = 0;
  
  RgKey nuclkey = "NuclearModel";
  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);  
}
//___________________________________________________________________________

