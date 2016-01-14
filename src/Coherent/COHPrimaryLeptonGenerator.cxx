//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the new Coherent package from its previous location  (EVGModules 
   package)

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TVector3.h>

#include "Conventions/Constants.h"
#include "Coherent/COHPrimaryLeptonGenerator.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
COHPrimaryLeptonGenerator::COHPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::COHPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
COHPrimaryLeptonGenerator::COHPrimaryLeptonGenerator(string config) :
PrimaryLeptonGenerator("genie::COHPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
COHPrimaryLeptonGenerator::~COHPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void COHPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton in COH events

  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  // Look-up selected kinematics
  double Q2 = interaction->Kine().Q2(true);
  double y  = interaction->Kine().y(true);

  // Auxiliary params
  double Ev  = init_state.ProbeE(kRfLab);
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  double El  = (1-y)*Ev;
  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)

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

  // Lepton 4-momentum in LAB
  TLorentzVector p4l(p3l,El);

  // Figure out the Final State Lepton PDG Code
  int pdgc = interaction->FSPrimLepton()->PdgCode();

  // Create a GHepParticle and add it to the event record
  this->AddToEventRecord(evrec, pdgc, p4l);

  // Set final state lepton polarization
  this->SetPolarization(evrec);
}
//___________________________________________________________________________
