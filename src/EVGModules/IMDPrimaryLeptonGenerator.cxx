//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - July 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TVector3.h>

#include "Conventions/Constants.h"
#include "EVGModules/IMDPrimaryLeptonGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
IMDPrimaryLeptonGenerator::IMDPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::IMDPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
IMDPrimaryLeptonGenerator::IMDPrimaryLeptonGenerator(string config):
PrimaryLeptonGenerator("genie::IMDPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
IMDPrimaryLeptonGenerator::~IMDPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void IMDPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton for IMD events

  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  // Get selected kinematics
  double y = interaction->Kine().y(true);
  assert(y>0 && y<1);

  // Final state primary lepton PDG code
  int pdgc = kPdgMuon;

  // Compute the neutrino and muon energy
  double Ev    = init_state.ProbeE(kRfLab); 
  double Emu   = (1-y)*Ev;

  // Compute the momentum transfer and scattering angle
  double Emu2  = TMath::Power(Emu,2);
  double me    = kElectronMass;
  double mmu2  = kMuonMass2;
  double pmu   = TMath::Sqrt(Emu2-mmu2);   
  
  assert(Emu2>=mmu2);

  double Q2    = 2*(Ev-Emu)*me;
  double costh = (Emu-0.5*(Q2+mmu2)/Ev)/pmu;
  double sinth = TMath::Sqrt( TMath::Max(0., 1-TMath::Power(costh,2.)) );

  //warn about overflow in costheta and ignore it if it is small or abort
  if( TMath::Abs(costh)>1 ) {
     LOG("LeptonicVertex", pWARN)
       << "El = " << Emu << ", Ev = " << Ev << ", cos(theta) = " << costh;
     if(TMath::Abs(costh)-1<0.3) costh = 1.0; //why?
  }
  assert(TMath::Abs(costh)<=1);

  // Compute the p components along and perpendicular the v direction 
  double plp = pmu * costh; // p(//)
  double plt = pmu * sinth; // p(-|)

  LOG("LeptonicVertex", pNOTICE)
        << "fsl: E = " << Emu << ", |p//| = " << plp << "[pT] = " << plt;

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

  // Lepton 4-momentum in the LAB
  TLorentzVector p4l(p3l,Emu);

  // Create a GHepParticle and add it to the event record
  this->AddToEventRecord(evrec, pdgc, p4l);

  // Set final state lepton polarization
  this->SetPolarization(evrec);
}
//___________________________________________________________________________
