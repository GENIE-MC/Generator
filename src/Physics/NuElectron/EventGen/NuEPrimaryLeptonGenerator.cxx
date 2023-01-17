//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

 Changes required to boost from electron rest frame to lab frame
 were installed by Brinden Carlson (Univ. of Florida)
*/
//____________________________________________________________________________

#include <TMath.h>
#include <TVector3.h>

#include "Framework/Conventions/Constants.h"
#include "Physics/NuElectron/EventGen/NuEPrimaryLeptonGenerator.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
NuEPrimaryLeptonGenerator::NuEPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::NuEPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
NuEPrimaryLeptonGenerator::NuEPrimaryLeptonGenerator(string config):
PrimaryLeptonGenerator("genie::NuEPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
NuEPrimaryLeptonGenerator::~NuEPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void NuEPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton for NuE events

  Interaction * interaction = evrec->Summary();

  // Boost vector for [LAB] <-> [Electron Rest Frame] transforms
  TVector3 beta = this->EleRestFrame2Lab(evrec); // Get boost of hit

  // Neutrino 4p
  TLorentzVector * p4v = evrec->Probe()->GetP4(); // v 4p @ LAB
  p4v->Boost(-1.*beta);                           // v 4p @ Electron rest frame

  // Get selected kinematics
  double y = interaction->Kine().y(true);
  assert(y>0 && y<1);

  // Final state primary lepton PDG code
  int pdgc = interaction->FSPrimLeptonPdg();
  assert(pdgc!=0);

  // Compute the neutrino and muon energy
  double Ev  = p4v->E();
  double El  = (1-y)*Ev;

  LOG("LeptonicVertex", pINFO)
       << "Ev = " << Ev << ", y = " << y << ", -> El = " << El;

  // Compute the momentum transfer and scattering angle
  double El2   = TMath::Power(El,2);
  double me    = kElectronMass;
  double ml    = interaction->FSPrimLepton()->Mass();
  double ml2   = TMath::Power(ml,2);
  double pl    = TMath::Sqrt(El2-ml2);

  assert(El2>=ml2);

  double Q2    = 2*(Ev-El)*me + me*me;
  double costh = (El-0.5*(Q2+ml2)/Ev)/pl;
  double sinth = TMath::Sqrt( TMath::Max(0., 1-TMath::Power(costh,2.)) );

  LOG("LeptonicVertex", pNOTICE)
       << "Q2 = " << Q2 << ", cos(theta) = " << costh;

  //warn about overflow in costheta and ignore it if it is small or abort
  if( TMath::Abs(costh)>1 ) {
     LOG("LeptonicVertex", pWARN)
       << "El = " << El << ", Ev = " << Ev << ", cos(theta) = " << costh;
     //if(TMath::Abs(costh)-1<0.3) costh = 1.0; //why?
  }
  assert(TMath::Abs(costh)<=1);

  // Compute the p components along and perpendicular the v direction
  double plp = pl * costh; // p(//)
  double plt = pl * sinth; // p(-|)

  LOG("LeptonicVertex", pNOTICE)
        << "fsl: E = " << El << ", |p//| = " << plp << "[pT] = " << plt;

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);

  // Take a unit vector along the neutrino direction
  TVector3 unit_nudir = evrec->Probe()->P4()->Vect().Unit();

  // Rotate lepton momentum vector from the reference frame (x'y'z') where
  // {z':(neutrino direction), z'x':(theta plane)} to the LAB
  TVector3 p3l(pltx,plty,plp);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in the Ele rest frame
  TLorentzVector p4l(p3l,El);
  p4l.Boost(beta); //Boost back to lab

  // Create a GHepParticle and add it to the event record
  this->AddToEventRecord(evrec, pdgc, p4l);

  // Set final state lepton polarization
  this->SetPolarization(evrec);
}
//___________________________________________________________________________
TVector3 NuEPrimaryLeptonGenerator::EleRestFrame2Lab(GHepRecord * evrec) const
{
// Velocity for an active Lorentz transform taking the final state primary
// lepton from the [electron rest frame] --> [LAB]

  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  const TLorentzVector & pele4 = init_state.Tgt().HitEleP4(); //[@LAB]
  TVector3 beta = pele4.BoostVector();

  return beta;
}
//___________________________________________________________________________
