//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 03, 2008 - CA
   First added in v2.7.1

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
#include "Physics/NucleonDecay/NucleonDecayPrimaryVtxGenerator.h"
#include "Physics/NucleonDecay/NucleonDecayUtils.h"
#include "Physics/NucleonDecay/NucleonDecayMode.h"
#include "Physics/NuclearState/NuclearModelI.h"

using namespace genie;

//____________________________________________________________________________
NucleonDecayPrimaryVtxGenerator::NucleonDecayPrimaryVtxGenerator() :
EventRecordVisitorI("genie::NucleonDecayPrimaryVtxGenerator")
{

}
//____________________________________________________________________________
NucleonDecayPrimaryVtxGenerator::NucleonDecayPrimaryVtxGenerator(
  string config) :
EventRecordVisitorI("genie::NucleonDecayPrimaryVtxGenerator",config)
{

}
//____________________________________________________________________________
NucleonDecayPrimaryVtxGenerator::~NucleonDecayPrimaryVtxGenerator() 
{

}
//____________________________________________________________________________
void NucleonDecayPrimaryVtxGenerator::ProcessEventRecord(
  GHepRecord * event) const
{
  Interaction * interaction = event->Summary();
  fCurrInitStatePdg = interaction->InitState().Tgt().Pdg();
  fCurrDecayMode = (NucleonDecayMode_t) interaction->ExclTag().DecayMode();
  fCurrDecayedNucleon = interaction->InitState().Tgt().HitNucPdg();

  LOG("NucleonDecay", pNOTICE)
    << "Simulating decay " << utils::nucleon_decay::AsString(fCurrDecayMode, fCurrDecayedNucleon)
    << " for an initial state with code: " << fCurrInitStatePdg;

  fNucleonIsBound = (pdg::IonPdgCodeToA(fCurrInitStatePdg) > 1);

  this->AddInitialState(event);
  this->GenerateDecayedNucleonPosition(event);
  this->GenerateFermiMomentum(event);
  this->GenerateDecayProducts(event);
}
//____________________________________________________________________________
void NucleonDecayPrimaryVtxGenerator::AddInitialState(
  GHepRecord * event) const
{
//
// Add initial state in the event record.
//
// If the decayed nucleon is one bound in a nucleus, the event record is
// initialized as follows:
//    id: 0, nucleus (kIStInitialState)
//    |     
//    |---> { id: 1, nucleon         (kIStDecayedState)
//          { id: 2, remnant nucleus (kIStStableFinalState)
//
// If the decayed nucleon is a free one, the event record is initialized as
// follows:
//    id: 0, nucleon (kIStInitialState)
//    |     
//    |---> id: 1, nucleon (kIStDecayedState)
//

  TLorentzVector v4(0,0,0,0);
  
  GHepStatus_t stis = kIStInitialState;
  GHepStatus_t stdc = kIStDecayedState;
  GHepStatus_t stfs = kIStStableFinalState;

  int ipdg = fCurrInitStatePdg;
  
  // Decayed nucleon is a bound one.
  if(fNucleonIsBound) 
  {
    // add initial nucleus
    double Mi  = PDGLibrary::Instance()->Find(ipdg)->Mass();
    TLorentzVector p4i(0,0,0,Mi);
    event->AddParticle(ipdg,stis,-1,-1,-1,-1, p4i, v4);
               
    // add decayed nucleon
    int dpdg = fCurrDecayedNucleon;
    double mn = PDGLibrary::Instance()->Find(dpdg)->Mass();
    TLorentzVector p4n(0,0,0,mn);  
    event->AddParticle(dpdg,stdc, 0,-1,-1,-1, p4n, v4);
     
    // add nuclear remnant
    int A = pdg::IonPdgCodeToA(ipdg);
    int Z = pdg::IonPdgCodeToZ(ipdg);
    A--;
    if(dpdg == kPdgProton) { Z--; }
    int rpdg = pdg::IonPdgCode(A, Z);
    double Mf  = PDGLibrary::Instance()->Find(rpdg)->Mass();
    TLorentzVector p4f(0,0,0,Mf);
    event->AddParticle(rpdg,stfs,0,-1,-1,-1, p4f, v4);
  }

  // Decayed nucleon is a free one
  else
  {
    // Initial state is either a neutron or a proton.
    // Convert the initial state PDG code from the ion convention (10LZZZAAAI)
    // to the usual code for neutrons or protons
    int ipdg_short = -1;
    if(ipdg == kPdgTgtFreeP) ipdg_short = kPdgProton;
    if(ipdg == kPdgTgtFreeN) ipdg_short = kPdgNeutron;

    // Decayed nucleon code 
    int dpdg = fCurrDecayedNucleon;

    if(dpdg != ipdg_short) {
       LOG("NucleonDecay", pWARN)
           << "Couldn't generate given decay (" 
           << utils::nucleon_decay::AsString(fCurrDecayMode, fCurrDecayedNucleon) << ")"
           << " for given initial state (PDG = " << ipdg_short << ")";
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Decay-mode / Initial-state mismatch!");
       exception.SwitchOnFastForward();
       throw exception;
    }
    // add initial nucleon
    double mn  = PDGLibrary::Instance()->Find(ipdg)->Mass();
    TLorentzVector p4i(0,0,0,mn);
    event->AddParticle(dpdg,stis,-1,-1,-1,-1, p4i, v4);
    // add decayed nucleon
    event->AddParticle(dpdg,stdc,0,-1,-1,-1, p4i, v4);
  }
}
//____________________________________________________________________________
void NucleonDecayPrimaryVtxGenerator::GenerateDecayedNucleonPosition(
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
            
  LOG("NucleonDecay", pINFO)
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
       LOG("NucleonDecay", pWARN)
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
       LOG("NucleonDecay", pERROR)
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

  GHepParticle * decayed_nucleon = event->Particle(1);
  assert(decayed_nucleon);
  decayed_nucleon->SetPosition(vtx);
}
//____________________________________________________________________________
void NucleonDecayPrimaryVtxGenerator::GenerateFermiMomentum(
  GHepRecord * event) const
{
  GHepParticle * initial_nucleus = event->Particle(0);
  int A = initial_nucleus->A();
  if(A <= 2) {
    return;
  }

  GHepParticle * decayed_nucleon = event->Particle(1);
  GHepParticle * remnant_nucleus = event->Particle(2);
  assert(decayed_nucleon);
  assert(remnant_nucleus);

  // generate a Fermi momentum & removal energy
  Target tgt(initial_nucleus->Pdg());
  tgt.SetHitNucPdg(decayed_nucleon->Pdg());
  fNuclModel->GenerateNucleon(tgt);
  TVector3 p3 = fNuclModel->Momentum3();
  double w    = fNuclModel->RemovalEnergy();

  LOG("FermiMover", pINFO) 
     << "Generated nucleon momentum: ("
     << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << "), "
     << "|p| = " << p3.Mag();
  LOG("NucleonDecay", pINFO) 
     << "Generated nucleon removal energy: w = " << w;
  
  double pF2 = p3.Mag2(); // (fermi momentum)^2

  double Mi  = PDGLibrary::Instance()->Find(initial_nucleus->Pdg())-> Mass(); // initial nucleus mass
  double Mf  = PDGLibrary::Instance()->Find(remnant_nucleus->Pdg())-> Mass(); // remnant nucleus mass

  double EN = Mi - TMath::Sqrt(pF2 + Mf*Mf);

  TLorentzVector p4nucleon(   p3.Px(),    p3.Py(),    p3.Pz(), EN);
  TLorentzVector p4remnant(-1*p3.Px(), -1*p3.Py(), -1*p3.Pz(), Mi-EN);

  decayed_nucleon->SetMomentum(p4nucleon);
  remnant_nucleus->SetMomentum(p4remnant);
}
//____________________________________________________________________________
void NucleonDecayPrimaryVtxGenerator::GenerateDecayProducts(
  GHepRecord * event) const
{
  LOG("NucleonDecay", pINFO) << "Generating decay...";

  PDGCodeList pdgv = utils::nucleon_decay::DecayProductList(fCurrDecayMode, fCurrDecayedNucleon);
  LOG("NucleonDecay", pINFO) << "Decay product IDs: " << pdgv;
  assert ( pdgv.size() >  1);

  LOG("NucleonDecay", pINFO) << "Performing a phase space decay...";

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

  LOG("NucleonDecay", pINFO)  
    << "Decaying N = " << pdgv.size() << " particles / total mass = " << sum;

  int decayed_nucleon_id = 1;
  GHepParticle * decayed_nucleon = event->Particle(decayed_nucleon_id);
  assert(decayed_nucleon);
  TLorentzVector * p4d = decayed_nucleon->GetP4();
  TLorentzVector * v4d = decayed_nucleon->GetX4();

  LOG("NucleonDecay", pINFO) 
    << "Decaying system p4 = " << utils::print::P4AsString(p4d);

  // Set the decay
  bool permitted = fPhaseSpaceGenerator.SetDecay(*p4d, pdgv.size(), mass);
  if(!permitted) {
     LOG("NucleonDecay", pERROR) 
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

  LOG("NucleonDecay", pNOTICE) 
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
       LOG("NucleonDecay", pWARN) 
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
        LOG("NucleonDecay", pWARN) 
           << "Decay weight = " << w << " > max decay weight = " << wmax;
     }
     double gw = wmax * rnd->RndHadro().Rndm();
     accept_decay = (gw<=w);

     LOG("NucleonDecay", pINFO) 
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
        utils::nucleon_decay::DecayProductStatus(fNucleonIsBound, pdgc);
     event->AddParticle(pdgc, ist, decayed_nucleon_id,-1,-1,-1, *p4fin, v4);
     idp++;
  }

  // Clean-up
  delete [] mass;
  delete p4d;
  delete v4d;
}
//____________________________________________________________________________
void NucleonDecayPrimaryVtxGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);   
  this->LoadConfig();
}
//___________________________________________________________________________
void NucleonDecayPrimaryVtxGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void NucleonDecayPrimaryVtxGenerator::LoadConfig(void)
{
//  AlgConfigPool * confp = AlgConfigPool::Instance();
//  const Registry * gc = confp->GlobalParameterList();
    
  fNuclModel = 0;
  
  RgKey nuclkey = "NuclearModel";
  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);  
}
//___________________________________________________________________________

