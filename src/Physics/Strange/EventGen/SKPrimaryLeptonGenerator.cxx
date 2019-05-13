//____________________________________________________________________________
/*

  Copyright (c) 2003-2019, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE
  
  Authors: Chris Marshall <marshall \at pas.rochester.edu>
           University of Rochester

           Martti Nirkko
           University of Berne
 
   For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/Strange/EventGen/SKPrimaryLeptonGenerator.h"

using namespace genie;

//___________________________________________________________________________
SKPrimaryLeptonGenerator::SKPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::SKPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
SKPrimaryLeptonGenerator::SKPrimaryLeptonGenerator(string config) :
PrimaryLeptonGenerator("genie::SKPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
SKPrimaryLeptonGenerator::~SKPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void SKPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{

  // no modification is required to the std implementation
  PrimaryLeptonGenerator::ProcessEventRecord(evrec);

  if(evrec->FinalStatePrimaryLepton()->IsOffMassShell()) {
    LOG("LeptonicVertex", pERROR)
               << "*** Selected kinematics lead to off mass shell lepton!";
     evrec->EventFlags()->SetBitNumber(kLeptoGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("E<m for final state lepton");
     exception.SwitchOnFastForward();
     throw exception;
  }
  //CalculatePrimaryLepton(evrec);
}
/*
//___________________________________________________________________________
void SKPrimaryLeptonGenerator::CalculatePrimaryLepton(GHepRecord * evrec) const
{
// This method generates the final state primary lepton in single-K events

  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  // Look-up selected kinematics
  double lep_t = interaction->Kine().GetKV(kKVSelTl);
  double lep_costheta = interaction->Kine().GetKV(kKVSelctl);

  // Auxiliary params
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);
  double El = lep_t + ml;
  double plep = TMath::Sqrt( El*El - ml2 );

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();

  LOG( "SKLepton", pDEBUG )
    << "lepton T = " << lep_t << " cos theta = " << lep_costheta << " random phi = " << phi;

  // Lepton 3vector w.r.t. neutrino direction
  TVector3 p3l(0,0,0);
  p3l.SetMagThetaPhi(plep, TMath::ACos(lep_costheta), phi);

  // Take a unit vector along the neutrino direction
  TVector3 unit_nudir = evrec->Probe()->P4()->Vect().Unit();

  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the LAB
  p3l.RotateUz(unit_nudir);

  LOG( "SKLepton", pDEBUG )
    << "lab frame lepton px = " << p3l.x() << " py = " << p3l.y() << " pz = " << p3l.z();

  // Lepton 4-momentum in LAB
  TLorentzVector p4l(p3l,El);

  // Figure out the Final State Lepton PDG Code
  int pdgc = interaction->FSPrimLepton()->PdgCode();

  // Create a GHepParticle and add it to the event record
  this->AddToEventRecord(evrec, pdgc, p4l);

  // Set final state lepton polarization
  this->SetPolarization(evrec);
}
*/
