//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TVector3.h>
#include <TLorentzVector.h>

#include "Conventions/Constants.h"
#include "EVGModules/PrimaryLeptonGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "Numerical/RandomGen.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
PrimaryLeptonGenerator::PrimaryLeptonGenerator() :
EventRecordVisitorI()
{

}
//___________________________________________________________________________
PrimaryLeptonGenerator::PrimaryLeptonGenerator(string name) :
EventRecordVisitorI(name)
{

}
//___________________________________________________________________________
PrimaryLeptonGenerator::PrimaryLeptonGenerator(string name, string config) :
EventRecordVisitorI(name, config)
{

}
//___________________________________________________________________________
PrimaryLeptonGenerator::~PrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void PrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton

  Interaction * interaction = evrec->GetInteraction();

  // Boost vector for [LAB] <-> [Nucleon Rest Frame] transforms
  TVector3 beta = this->NucRestFrame2Lab(evrec);

  // Neutrino 4p
  TLorentzVector * p4v = evrec->Probe()->GetP4(); // v 4p @ LAB
  p4v->Boost(-1.*beta);                           // v 4p @ Nucleon rest frame

  // Look-up selected kinematics & other needed kinematical params
  double Q2  = interaction->GetKinematics().Q2(true);
  double y   = interaction->GetKinematics().y(true);
  double Ev  = p4v->E(); 
  double ml  = interaction->GetFSPrimaryLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  LOG("LeptonicVertex", pNOTICE)
             << "Ev = " << Ev << ", Q2 = " << Q2 << ", y = " << y;

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  double El  = (1-y)*Ev;
  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)

  LOG("LeptonicVertex", pNOTICE)
        << "fsl: E = " << El << ", |p//| = " << plp << ", [pT] = " << plt;

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);

  // Take a unit vector along the neutrino direction @ the nucleon rest frame
  TVector3 unit_nudir = p4v->Vect().Unit(); 

  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the nucleon rest frame
  TVector3 p3l(pltx,plty,plp);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in the nucleon rest frame
  TLorentzVector p4l(p3l,El);

  LOG("LeptonicVertex", pNOTICE)
       << "fsl @ NRF: " << utils::print::P4AsString(&p4l);

  // Boost final state primary lepton to the lab frame
  p4l.Boost(beta); // active Lorentz transform

  LOG("LeptonicVertex", pNOTICE)
       << "fsl @ LAB: " << utils::print::P4AsString(&p4l);

  // Figure out the Final State Lepton PDG Code
  int pdgc = interaction->GetFSPrimaryLepton()->PdgCode();

  // Create a GHepParticle and add it to the event record
  this->AddToEventRecord(evrec, pdgc, p4l);

  // Set final state lepton polarization
  this->SetPolarization(evrec);

  delete p4v;
}
//___________________________________________________________________________
TVector3 PrimaryLeptonGenerator::NucRestFrame2Lab(GHepRecord * evrec) const
{
// Velocity for an active Lorentz transform taking the final state primary
// lepton from the [nucleon rest frame] --> [LAB]

  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  TLorentzVector * pnuc4 = init_state.GetTarget().StruckNucleonP4(); //[LAB]

  double bx = pnuc4->Px() / pnuc4->Energy();
  double by = pnuc4->Py() / pnuc4->Energy();
  double bz = pnuc4->Pz() / pnuc4->Energy();

  TVector3 beta(bx,by,bz);
  return beta;
}
//___________________________________________________________________________
void PrimaryLeptonGenerator::AddToEventRecord(
              GHepRecord * evrec, int pdgc, const TLorentzVector & p4) const
{
// Adds the final state primary lepton GHepParticle to the event record.
// To be called by all concrete PrimaryLeptonGenerators before exiting.

  int mom = evrec->ProbePosition();
  TLorentzVector vdummy(0,0,0,0); // position 4-vector

  evrec->AddParticle(pdgc, kIStStableFinalState, mom,-1,-1,-1, p4, vdummy);

  // update the interaction summary
  evrec->GetInteraction()->GetKinematicsPtr()->SetFSLeptonP4(p4);
}
//___________________________________________________________________________
void PrimaryLeptonGenerator::SetPolarization(GHepRecord * ev) const
{
// Set the final state lepton polarization. A mass-less lepton would be fully
// polarized. This would be exact for neutrinos and a very good approximation
// for electrons for the energies this generator is going to be used. This is
// not the case for muons and, mainly, for taus. I need to refine this later.
// How? See Kuzmin, Lyubushkin and Naumov, hep-ph/0312107

  // get the final state primary lepton
  GHepParticle * fsl = ev->FinalStatePrimaryLepton();
  if(!fsl) {
     LOG("LeptonicVertex", pERROR)
                    << "Final state lepton not set yet! \n" << *ev;
     return;
  }

  // get (px,py,pz) @ LAB
  TVector3 plab(fsl->Px(), fsl->Py(), fsl->Pz());

  // in the limit m/E->0: leptons are left-handed and their anti-particles
  // are right-handed
  int pdgc = fsl->PdgCode();
  if(pdg::IsNeutrino(pdgc) || pdg::IsElectron(pdgc) ||
                    pdg::IsMuon(pdgc) || pdg::IsTau(pdgc) ) {
    plab *= -1; // left-handed
  }

  LOG("LeptonicVertex", pINFO) 
            << "Setting polarization angles for particle: " << fsl->Name();

  fsl->SetPolarization(plab);

  if(fsl->PolzIsSet()) {
     LOG("LeptonicVertex", pINFO) 
          << "Polarization (rad): Polar = "  << fsl->PolzPolarAngle() 
                           << ", Azimuthal = " << fsl->PolzAzimuthAngle();
  }
}
//___________________________________________________________________________

