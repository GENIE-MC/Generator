//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 12, 2013 - CA (code from Rosen Matev)
   Handle the IMD annihilation channel.

*/
//____________________________________________________________________________

#include <TVector3.h>
#include <TLorentzVector.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Physics/Common/PrimaryLeptonGenerator.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Numerical/MathUtils.h"

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

  Interaction * interaction = evrec->Summary();

  // Boost vector for [LAB] <-> [Nucleon Rest Frame] transforms
  TVector3 beta = this->NucRestFrame2Lab(evrec);

  // Neutrino 4p
  TLorentzVector * p4v = evrec->Probe()->GetP4(); // v 4p @ LAB
  p4v->Boost(-1.*beta);                           // v 4p @ Nucleon rest frame

  // Look-up selected kinematics & other needed kinematical params
  double Q2  = interaction->Kine().Q2(true);
  double y   = interaction->Kine().y(true);
  double Ev  = p4v->E(); 
  double ml  = interaction->FSPrimLepton()->Mass();
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
  int pdgc = interaction->FSPrimLepton()->PdgCode();

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

  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  const TLorentzVector & pnuc4 = init_state.Tgt().HitNucP4(); //[@LAB]
  TVector3 beta = pnuc4.BoostVector();
    
  return beta;
}
//___________________________________________________________________________
void PrimaryLeptonGenerator::AddToEventRecord(
              GHepRecord * evrec, int pdgc, const TLorentzVector & p4) const
{
// Adds the final state primary lepton GHepParticle to the event record.
// To be called by all concrete PrimaryLeptonGenerators before exiting.

  Interaction * interaction = evrec->Summary();
    
  GHepParticle * mom  = evrec->Probe();
  int            imom = evrec->ProbePosition();

  const TLorentzVector & vtx = *(mom->X4());

  TLorentzVector x4l(vtx);  // position 4-vector
  TLorentzVector p4l(p4);   // momentum 4-vector

  GHepParticle * nucltgt = evrec->TargetNucleus();

  bool is_ve = interaction->ProcInfo().IsInverseMuDecay() || 
               interaction->ProcInfo().IsIMDAnnihilation() || 
               interaction->ProcInfo().IsNuElectronElastic();

  bool can_correct = fApplyCoulombCorrection && nucltgt!=0 && !is_ve;
  if(can_correct) {
    LOG("LeptonicVertex", pINFO)  
        << "Correcting f/s lepton energy for Coulomb effects";

    double m = interaction->FSPrimLepton()->Mass();
    double Z  = nucltgt->Z();
    double A  = nucltgt->A();

    //  charge radius of nucleus in GeV^-1
    double Rc = (1.1*TMath::Power(A,1./3.) + 0.86*TMath::Power(A,-1./3.))/0.197; 

    // shift of lepton energy in homogenous sphere with radius Rc
    double Vo = 3*kAem*Z/(2*Rc);
    Vo *= 0.75; // as suggested in R.Gran's note
    
    double Elo = p4l.Energy();
    double e   = TMath::Min(Vo, Elo-m);
    double El  = TMath::Max(0., Elo-e);

    LOG("LeptonicVertex", pINFO) 
      << "Lepton energy subtraction: E = " << Elo << " --> " << El;

    double pmag_old = p4l.P();
    double pmag_new = TMath::Sqrt(utils::math::NonNegative(El*El-m*m));
    double scale    = pmag_new / pmag_old;
    LOG("LeptonicVertex", pDEBUG) 
         << "|pnew| = " << pmag_new << ", |pold| = " << pmag_old
         << ", scale = " << scale;

    double pxl = scale * p4l.Px();
    double pyl = scale * p4l.Py();
    double pzl = scale * p4l.Pz();

    p4l.SetPxPyPzE(pxl,pyl,pzl,El);

    TLorentzVector p4c = p4 - p4l;
    TLorentzVector x4dummy(0,0,0,0);;

    evrec->AddParticle(
        kPdgCoulobtron, kIStStableFinalState, -1,-1,-1,-1, p4c, x4dummy);
  }

  evrec->AddParticle(pdgc, kIStStableFinalState, imom,-1,-1,-1, p4l, x4l);

  // update the interaction summary
  evrec->Summary()->KinePtr()->SetFSLeptonP4(p4l);
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
  int pdgc = fsl->Pdg();
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
void PrimaryLeptonGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PrimaryLeptonGenerator::Configure(string config)
{    
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PrimaryLeptonGenerator::LoadConfig(void)
{
  GetParam( "ApplyCoulombCorrection", fApplyCoulombCorrection ) ;

}
//___________________________________________________________________________

