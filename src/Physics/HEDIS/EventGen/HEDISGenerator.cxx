//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/HEDIS/EventGen/HEDISGenerator.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Conventions/Constants.h"

using namespace genie;
using namespace genie::utils::math;

//___________________________________________________________________________
HEDISGenerator::HEDISGenerator() :
HadronicSystemGenerator("genie::HEDISGenerator")
{
  this->Initialize();
}
//___________________________________________________________________________
HEDISGenerator::HEDISGenerator(string config) :
HadronicSystemGenerator("genie::HEDISGenerator", config)
{
  this->Initialize();
}
//___________________________________________________________________________
HEDISGenerator::~HEDISGenerator()
{

}
//____________________________________________________________________________
void HEDISGenerator::Initialize(void) const
{

}
//___________________________________________________________________________
void HEDISGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  //-- Add the target remnant
  this->AddTargetNucleusRemnant(evrec);
  GHepParticle * target = evrec -> TargetNucleus();
  if(target) evrec->Particle(evrec->RemnantNucleusPosition())->SetStatus(kIStFinalStateNuclearRemnant);

  //-- Add the primary lepton
  this->AddPrimaryLepton(evrec);

  //-- Run the hadronization model and get the fragmentation products
  fHadronizationModel->ProcessEventRecord(evrec);

}
//___________________________________________________________________________
void HEDISGenerator::AddPrimaryLepton(GHepRecord * evrec) const
{

  Interaction * interaction = evrec->Summary();

  // Neutrino 4p
  LongLorentzVector p4v( * evrec->Probe()->P4() );
  LOG("HEDISGenerator", pINFO) << "NEUTRINO @ LAB' =>  E = " << p4v.E() << " //  m = " << p4v.M() << " // p = " << p4v.P();
  LOG("HEDISGenerator", pINFO) << "                  dir = " << p4v.Dx() << " , "  << p4v.Dy() << " , "  << p4v.Dz();

  // Look-up selected kinematics & other needed kinematical params
  long double Q2  = interaction->Kine().Q2(true);
  long double y   = interaction->Kine().y(true);
  long double Ev  = p4v.E();
  long double ml  = interaction->FSPrimLepton()->Mass();
  long double ml2 = powl(ml,2);

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction
  long double El  = (1-y)*Ev;
  long double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  long double plt = sqrtl(fmaxl(0.,El*El-plp*plp-ml2)); // p(-|)
  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  long double phi  = 2 * constants::kPi * rnd->RndLep().Rndm();
  long double pltx = plt * cosl(phi);
  long double plty = plt * sinl(phi);

  // Lepton 4-momentum in the LAB frame
  LongLorentzVector p4llong( pltx, plty, plp, El );
  p4llong.Rotate(p4v);
  LOG("HEDISGenerator", pINFO) << "LEPTON     @ LAB' =>  E = " << p4llong.E() << " //  m = " << p4llong.M() << " // p = " << p4llong.P();
  LOG("HEDISGenerator", pINFO) << "                    dir = " << p4llong.Dx() << " , "  << p4llong.Dy() << " , "  << p4llong.Dz();

  // Translate from long double to double
  TLorentzVector p4l( (double)p4llong.Px(), (double)p4llong.Py(), (double)p4llong.Pz(), (double)p4llong.E() );

  // Add lepton to EventRecord
  int pdgl = interaction->FSPrimLepton()->PdgCode();
  evrec->AddParticle(pdgl, kIStStableFinalState, evrec->ProbePosition(),-1,-1,-1, p4l, *(evrec->Probe()->X4()));
  evrec->Summary()->KinePtr()->SetFSLeptonP4(p4l);

}
//___________________________________________________________________________
void HEDISGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISGenerator::LoadConfig(void)
{
  fHadronizationModel = 0;

  //-- Get the requested hadronization model
  fHadronizationModel =
     dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("Hadronizer"));
  assert(fHadronizationModel);

}
