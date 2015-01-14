//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - September 26, 2005

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
#include "AtharSingleKaon/ASKPrimaryLeptonGenerator.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Base/XSecAlgorithmI.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/RunningThreadInfo.h"


using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
ASKPrimaryLeptonGenerator::ASKPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::ASKPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
ASKPrimaryLeptonGenerator::ASKPrimaryLeptonGenerator(string config) :
PrimaryLeptonGenerator("genie::ASKPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
ASKPrimaryLeptonGenerator::~ASKPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void ASKPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
//-- Access cross section algorithm for running thread
  //RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  //const EventGeneratorI * evg = rtinfo->RunningThread();
  //const XSecAlgorithmI *fXSecModel = evg->CrossSectionAlg();
  CalculatePrimaryLepton(evrec);
}
//___________________________________________________________________________
void ASKPrimaryLeptonGenerator::CalculatePrimaryLepton(GHepRecord * evrec) const
{
// This method generates the final state primary lepton in ASK events

  Interaction * interaction = evrec->Summary();
//  const InitialState & init_state = interaction->InitState();

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

  LOG( "ASKLepton", pDEBUG )
    << "lepton T = " << lep_t << " cos theta = " << lep_costheta << " random phi = " << phi;

  // Lepton 3vector w.r.t. neutrino direction
  TVector3 p3l(0,0,0);
  p3l.SetMagThetaPhi(plep, TMath::ACos(lep_costheta), phi);

  // Take a unit vector along the neutrino direction
  TVector3 unit_nudir = evrec->Probe()->P4()->Vect().Unit();

  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the LAB
  p3l.RotateUz(unit_nudir);

  LOG( "ASKLepton", pDEBUG )
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
