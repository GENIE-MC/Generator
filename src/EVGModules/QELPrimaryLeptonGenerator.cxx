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

#include "EVGModules/QELPrimaryLeptonGenerator.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"

using namespace genie;

//___________________________________________________________________________
QELPrimaryLeptonGenerator::QELPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::QELPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
QELPrimaryLeptonGenerator::QELPrimaryLeptonGenerator(string config):
PrimaryLeptonGenerator("genie::QELPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
QELPrimaryLeptonGenerator::~QELPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void QELPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton

  //-- Get the interaction & initial state objects
  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  //-- Figure out the Final State Lepton PDG Code
  int pdgc = interaction->GetFSPrimaryLepton()->PdgCode();

  //-- Use selected kinematics
  interaction->GetKinematicsPtr()->UseSelectedKinematics();

  //-- QEL Kinematics: Compute the lepton energy and the scattering
  //   angle with respect to the incoming neutrino

  //auxiliary params:
  double Ev   = init_state.GetProbeE(kRfStruckNucAtRest);
  double Q2   = interaction->GetKinematics().Q2();
  double M    = init_state.GetTarget().StruckNucleonP4()->M(); // can be off m/shell
  double ml   = interaction->GetFSPrimaryLepton()->Mass();
  double M2   = M*M;
  double ml2  = ml*ml;
  double s    = M2 + 2*M*Ev;
  double W2   = M2; // QEL: W ~= M
  double tmp  = s+ml2-W2;
  double tmp2 = tmp*tmp;

  //Compute outgoing lepton scat. angle with respect to the incoming v

  double cThSc = (tmp - (ml2+Q2)*2*s/(s-M2)) / (TMath::Sqrt(tmp2-4*s*ml2));
  assert( TMath::Abs(cThSc) <= 1 );

  //Compute outgoing lepton energy

  double El = Ev / ( 1 + (Ev/M) * (1-cThSc) );

  //-- Rotate its 4-momentum to the nucleon rest frame
  //   unit' = R(Theta0,Phi0) * R(ThetaSc,PhiSc) * R^-1(Theta0,Phi0) * unit

  TLorentzVector * p4l = P4InNucRestFrame(evrec, cThSc, El);

  //-- Boost it to the lab frame
  TVector3 * beta = NucRestFrame2Lab(evrec);
  p4l->Boost(*beta); // active Lorentz transform
  delete beta;

  //-- Create a GHepParticle and add it to the event record
  //   (use the insertion method at the base PrimaryLeptonGenerator visitor)

  this->AddToEventRecord(evrec, pdgc, p4l);

  delete p4l;

  //-- Set final state lepton polarization
  this->SetPolarization(evrec);

  //-- Reset running kinematics
  interaction->GetKinematicsPtr()->ClearRunningValues();
}
//___________________________________________________________________________
