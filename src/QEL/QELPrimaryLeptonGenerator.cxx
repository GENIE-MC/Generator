//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Moved into the new QEL package from its previous location (EVGModules)
 @ Mar 18, 2016 - JJ (SD)
   SetRunningLepton() is used by NievesQELCCPXSec::XSec() to generate a
   lepton before calculating the cross section. When Q2 is selected,
   QELKinematicsGenerator will store the generated lepton in the selected
   kinematics variables. If there is a stored selected lepton, this class
   adds that selected lepton to the event record. Otherwise the parent
   ProcessEventRecord() method is called.
*/
//____________________________________________________________________________

#include <TVector3.h>
#include <TLorentzVector.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Conventions/KineVar.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Interaction/Interaction.h"
#include "Interaction/Kinematics.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "QEL/QELPrimaryLeptonGenerator.h"
#include "Utils/PrintUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

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
  //PrimaryLeptonGenerator::ProcessEventRecord(evrec);
  Interaction * interaction = evrec->Summary();
  Kinematics * kine = interaction->KinePtr();

  // Check that a lepton was stored. If not, call the parent method.
  if(!(kine->KVSet(kKVSelTl) && kine->KVSet(kKVSelctl) 
       && kine->KVSet(kKVSelphikq))){
    LOG("LeptonicVertex",pINFO) << "Calling parent method to generate lepton";
    PrimaryLeptonGenerator::ProcessEventRecord(evrec);
    return;
  }
  LOG("LeptonicVertex",pINFO) << "Loading stored lepton";
  // Get the final state primary lepton energy and momentum components
  // along and perpendicular to the neutrino direction 
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  // Get the components stored by the QEL kinematics generator
  double Tl = kine->GetKV(kKVSelTl); // Used to store momentum magnitude pl
  double ctl = kine->GetKV(kKVSelctl);
  double phi = kine->GetKV(kKVSelphikq);

  double El = TMath::Sqrt(Tl*Tl + ml2);
  double stl = TMath::Sqrt(1-ctl*ctl);
  double plp = Tl * ctl;
  double pltx = Tl * stl * TMath::Cos(phi);
  double plty = Tl * stl * TMath::Sin(phi);

  TLorentzVector p4l(pltx,plty,plp,El);

  LOG("LeptonicVertex", pNOTICE)
       << "fsl @ LAB: " << utils::print::P4AsString(&p4l);

  // Figure out the Final State Lepton PDG Code
  int pdgc = interaction->FSPrimLepton()->PdgCode();

  // Create a GHepParticle and add it to the event record
  this->AddToEventRecord(evrec, pdgc, p4l);

  // Set final state lepton polarization
  this->SetPolarization(evrec);
}
//___________________________________________________________________________
// Generate a lepton using PrimaryLeptonGenerator::ProcessEventRecord()
// and store it by using the KineVar construct with
// the Kinematics class (variables kKVTl, kKVctl, kKVphikq).
// Stored lepton is in the LAB FRAME.
void QELPrimaryLeptonGenerator::SetRunningLepton(GHepRecord * evrec) const{
  Interaction * interaction = evrec->Summary();

  // Velocity for an active Lorentz transform taking the final state primary
  // lepton from the [nucleon rest frame] --> [LAB]
  const InitialState & init_state = interaction->InitState();
  const TLorentzVector & pnuc4 = init_state.Tgt().HitNucP4(); //[@LAB]
  TVector3 beta = pnuc4.BoostVector();


  // Neutrino 4p
  TLorentzVector * p4v = evrec->Probe()->GetP4(); // v 4p @ LAB
  p4v->Boost(-1.*beta);                           // v 4p @ Nucleon rest frame

  // Look-up selected kinematics & other needed kinematical params
  double Q2  = interaction->Kine().Q2(false);

  // get neutrino energy at struck nucleon rest frame and the
  // struck nucleon mass (can be off the mass shell)
  double E  = init_state.ProbeE(kRfHitNucRest);
  double M = init_state.Tgt().HitNucP4().M();

  const XclsTag & xcls = interaction->ExclTag();
  int rpdgc = 0;
  if(xcls.IsCharmEvent()) { rpdgc = xcls.CharmHadronPdg();           }
  else                    { rpdgc = interaction->RecoilNucleonPdg(); }
  assert(rpdgc);
  double W = PDGLibrary::Instance()->Find(rpdgc)->Mass();
  // (W,Q2) -> (x,y)
  double x=0, y=0;
  utils::kinematics::WQ2toXY(E,M,W,Q2,x,y);
  double Ev  = p4v->E(); 
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  LOG("LeptonicVertex", pINFO)
    << "Ev = " << Ev << ", Q2 = " << Q2 << ", y = " << y;

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  double El = (1-y)*Ev;
  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)

  LOG("LeptonicVertex", pINFO)
    << "trying fsl: E = " << El << ", |p//| = " << plp << ", [pT] = " << plt;

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

  LOG("LeptonicVertex", pINFO)
    << "trying fsl @ NRF: " << utils::print::P4AsString(&p4l);

  // Boost final state primary lepton to the lab frame
  p4l.Boost(beta); // active Lorentz transform

  LOG("LeptonicVertex", pINFO)
    << "trying fsl @ LAB: " << utils::print::P4AsString(&p4l);

  // Store the outgoing lepton components in the lab frame
  // Assume El^2 = pl^2 + ml^2
  double Tl = p4l.Vect().Mag(); // used to store 3-vector magnitude pl
  LOG("LeptonicVertex",pDEBUG) << "ml = " << ml << ", Tl = " << Tl;
  double ctl = p4l.CosTheta(); // cos(theta) in the lab frame
  phi = p4l.Phi(); // get phi in the lab fram
  interaction->KinePtr()->SetKV(kKVTl,Tl);
  interaction->KinePtr()->SetKV(kKVctl,ctl);
  interaction->KinePtr()->SetKV(kKVphikq,phi);

  delete p4v;
  return;

  /*LOG("LeptonicVertex",pFATAL) << "TESTING: evrec Q2 = " 
			       << evrec->Summary()->KinePtr()->Q2();
  // Copy the event record to generate a temporary lepton
  GHepRecord * temp_evrec = new GHepRecord(*evrec);
  Interaction * temp_in = new Interaction(*evrec->Summary());
  temp_evrec->AttachSummary(temp_in);

  // The hadronic inv. mass is equal to the recoil nucleon on-shell mass.
  // For QEL/Charm events it is set to be equal to the on-shell mass of
  // the generated charm baryon (Lamda_c+, Sigma_c+ or Sigma_c++)
  //
  const XclsTag & xcls = temp_in->ExclTag();
  int rpdgc = 0;
  if(xcls.IsCharmEvent()) { rpdgc = xcls.CharmHadronPdg();           }
  else                    { rpdgc = temp_in->RecoilNucleonPdg(); }
  assert(rpdgc);
  double gW = PDGLibrary::Instance()->Find(rpdgc)->Mass();

  const InitialState & init_state = temp_in->InitState();
  double E  = init_state.ProbeE(kRfHitNucRest);
  double M = init_state.Tgt().HitNucP4().M();
  
  double gQ2 = temp_evrec->Summary()->KinePtr()->Q2(false);
  // (W,Q2) -> (x,y)
  double gx=0, gy=0;
  utils::kinematics::WQ2toXY(E,M,gW,gQ2,gx,gy);

  temp_in->KinePtr()->SetQ2(gQ2,true);
  temp_in->KinePtr()->SetW(gW,true);
  temp_in->KinePtr()->Sety(gy,true);
  LOG("LeptonicVertex",pFATAL) << "TESTING: temp_evrec Q2 = " 
			       << temp_evrec->Summary()->KinePtr()->Q2(true);
  PrimaryLeptonGenerator::ProcessEventRecord(temp_evrec);

  // Store the outgoing lepton components in the lab frame
  // Assume El^2 = pl^2 + ml^2
  TLorentzVector p4l = temp_evrec->Summary()->KinePtr()->FSLeptonP4();
  double Tl = p4l.Vect().Mag(); // used to store 3-vector magnitude pl
  double ctl = p4l.CosTheta(); // cos(theta) in the lab frame
  double phi = p4l.Phi(); // get phi in the lab frame
  LOG("LeptonicVertex",pFATAL) << "TESTING: Tl = " << Tl << ", ctl = " << ctl
			       << ", phi = " << phi;
  Interaction * interaction = evrec->Summary();
  interaction->KinePtr()->SetKV(kKVTl,Tl);
  interaction->KinePtr()->SetKV(kKVctl,ctl);
  interaction->KinePtr()->SetKV(kKVphikq,phi);

  delete temp_evrec;
  return;*/
}
//___________________________________________________________________________
