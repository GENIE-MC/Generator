//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - September 26, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "EVGModules/COHPrimaryLeptonGenerator.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

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
// This method generates the final state primary lepton

  //-- Get the attached interaction & get initial state & kinematics objects
  Interaction * interaction = evrec->GetInteraction();

  const InitialState & init_state = interaction->GetInitialState();
  const Kinematics &   kinematics = interaction->GetKinematics();

  //-- Use selected kinematics
  interaction->GetKinematicsPtr()->UseSelectedKinematics();

  //-- Coherent Scattering Kinematics: Compute the lepton energy and the
  //   scattering angle with respect to the incoming neutrino

  //auxiliary vars
  TParticlePDG * fsl = interaction->GetFSPrimaryLepton();
  int    pdgc = fsl->PdgCode();
  double Ev   = init_state.GetProbeE(kRfLab);
  double x    = kinematics.x();
  double y    = kinematics.y();
  double ml   = fsl->Mass();
  double ml2  = ml*ml;
  double Q2   = 2*x*y*kNucleonMass*Ev;

  //Compute outgoing lepton energy & momentum
  double El = (1-y)*Ev;
  double pl = TMath::Sqrt( TMath::Max(0., El*El-ml2) );

  //Compute outgoing lepton scat. angle with respect to the incoming v
  double costheta = (El - 0.5*(Q2+ml2)/Ev) / pl; 

  LOG("LeptonicVertex", pDEBUG) << "cos(neutrino, fsl) = " << costheta;

  if(TMath::Abs(costheta) >= 1) {
     LOG("LeptonicVertex", pWARN) 
                  << "|cos(neutrino, fsl)| = " << costheta << " >= 1";
     costheta = TMath::Min(costheta,  1.);
     costheta = TMath::Max(costheta, -1.);
  }
  
  //-- Get the neutrino 4-p in LAB
  GHepParticle * neutrino = evrec->Probe();
  assert(neutrino);
  TLorentzVector * p4nu = neutrino->GetP4();

  //-- Rotate its 4-momentum to the LAB
  //   unit' = R(Theta0,Phi0) * R(ThetaSc,PhiSc) * R^-1(Theta0,Phi0) * unit
  TLorentzVector * p4l = this->Rotate4P(p4nu, pdgc, costheta, El);

  //-- Create a GHepParticle and add it to the event record
  //   (use the insertion method at the base PrimaryLeptonGenerator visitor)
  this->AddToEventRecord(evrec, pdgc, p4l);

  delete p4l;
  delete p4nu;

  // set final state lepton polarization
  this->SetPolarization(evrec);
}
//___________________________________________________________________________
