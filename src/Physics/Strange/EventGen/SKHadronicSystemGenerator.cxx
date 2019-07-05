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

#include <cstdlib>

#include <TVector3.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen//RunningThreadInfo.h"
#include "Physics/Strange/EventGen/SKHadronicSystemGenerator.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
SKHadronicSystemGenerator::SKHadronicSystemGenerator() :
HadronicSystemGenerator("genie::SKHadronicSystemGenerator")
{

}
//___________________________________________________________________________
SKHadronicSystemGenerator::SKHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::SKHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
SKHadronicSystemGenerator::~SKHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void SKHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Access cross section algorithm for running thread
  //RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  //const EventGeneratorI * evg = rtinfo->RunningThread();
  //const XSecAlgorithmI *fXSecModel = evg->CrossSectionAlg();
  CalculateHadronicSystem_AtharSingleKaon(evrec);
}
//___________________________________________________________________________
void SKHadronicSystemGenerator::CalculateHadronicSystem_AtharSingleKaon(GHepRecord * evrec) const
{
//
// This method generates the final state hadronic system (kaon + nucleus) 
//

  Interaction * interaction = evrec->Summary();
  Kinematics * kinematics = interaction->KinePtr();
  const XclsTag & xcls_tag  = interaction->ExclTag();

  //-- Access neutrino, initial nucleus and final state prim. lepton entries
  GHepParticle * nu  = evrec->Probe();
  GHepParticle * Ni  = evrec->HitNucleon();
  GHepParticle * fsl = evrec->FinalStatePrimaryLepton();
  assert(nu);
  assert(Ni);
  assert(fsl);

  const TLorentzVector vtx   = *(nu->X4());
  const TLorentzVector p4nu_lab  = *(nu ->P4());
  const TLorentzVector p4fsl_lab = *(fsl->P4());

  // these will get boosted to the nucleon rest frame
  TLorentzVector p4nu = p4nu_lab;
  TLorentzVector p4fsl = p4fsl_lab;

  // Transform the neutrino and final-state lepton to the struck nucleon rest frame
  const TLorentzVector pnuc4 = interaction->InitState().Tgt().HitNucP4(); // 4-momentum of struck nucleon in lab frame
  TVector3 beta = pnuc4.BoostVector();
  p4nu.Boost(-1.*beta);   
  p4fsl.Boost(-1.*beta);

  LOG( "SKHadron", pDEBUG ) << "\nStruck nucleon p = (" << pnuc4.X() << ", " << pnuc4.Y() << ", " << pnuc4.Z() << ")";
  LOG( "SKHadron", pDEBUG ) << "\nLab frame neutrino E = " << p4nu_lab.E() << " lepton " << p4fsl_lab.E() << " rest frame " << p4nu.E() << " lepton " << p4fsl.E();

  //-- Determine the pdg code of the final state nucleon
  int nuc_pdgc = (xcls_tag.NProtons()) ? kPdgProton : kPdgNeutron; // there's only ever one nucleon
  int kaon_pdgc = xcls_tag.StrangeHadronPdg();

  const Target & tgt = interaction->InitState().Tgt();
  GHepStatus_t ist = (tgt.IsNucleus()) ?
                          kIStHadronInTheNucleus : kIStStableFinalState;

  //-- basic kinematic inputs
  double Mf    = (xcls_tag.NProtons()) ? kProtonMass : kNeutronMass; // there's only ever one nucleon
  double M     = pnuc4.M();  // Mass of the struck nucleon
  double mk    = PDGLibrary::Instance()->Find(kaon_pdgc)->Mass(); // K+ and K0 mass are slightly different
  double mk2   = TMath::Power(mk,2);

  //-- specific kinematic quantities
  double kaon_T = kinematics->GetKV(kKVSelTk);
  double kaon_E = kaon_T + mk;
  double pk = sqrt( kaon_E*kaon_E - mk2 );

  TLorentzVector q = p4nu - p4fsl; // nucleon rest frame

  TVector3 qvec = q.Vect(); // nucleon rest frame
  double q3 = qvec.Mag(); // in q frame (0,0,q3)

  // Equation 17 of notes from M. Rafi Alam dated 6 November 2013
  double eN = q.E() + M - kaon_E; // FS nucleon total energy
  double cos_thetaKq = (q3*q3 + pk*pk + Mf*Mf - eN*eN)/(2*q3*pk);

  LOG( "SKHadron", pDEBUG ) << 
    "Cosine theta_kq = " << cos_thetaKq << "\n" <<
    "q.E = " << q.E() << " M = " << M << " kaon E " << kaon_E << " q3 = " << q3 << " pk = " << pk;

  if(cos_thetaKq > 1.0) {
     LOG("SKHadron", pWARN) << "Invalid selected kinematics; Attempt regenerating";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Invalid selected kinematics");
     exception.SwitchOnStepBack();
     exception.SetReturnStep(0);
     throw exception;
  }

//  this can be slightly larger than 1 due to numerical precision issues -- don't let it be
//  if( cos_thetaKq > 1.0 ) {
//    LOG( "SKHadron", pWARN ) << 
//      "Cosine theta_kq = " << cos_thetaKq << ", setting to 1.0\n" <<
//      "q.E = " << q.E() << " M = " << M << " kaon E " << kaon_E << " q3 = " << q3 << " pk = " << pk;
//    cos_thetaKq = 1.0;
//  }


  // Get phi for the k-q plane relative to nu-l plane
  double phi_kq = kinematics->GetKV(kKVSelphikq);

  TVector3 kaon;
  kaon.SetMagThetaPhi( pk, TMath::ACos(cos_thetaKq), phi_kq );

  TVector3 nucleon( -kaon.X(), -kaon.Y(), q3 - kaon.Z() ); // force 3-momentum conservation in q frame

  // Transform hadron momenta from z axis along q to lab frame
  kaon.RotateUz(qvec.Unit());
  nucleon.RotateUz(qvec.Unit());

  LOG( "SKHadron", pDEBUG ) <<
    "\nKaon (x,y,z) in nuc rest frame: (" << kaon.X() << ", " << kaon.Y() << ", " << kaon.Z() << ")" <<
    "\nNucleon:                        (" << nucleon.X() << ", " << nucleon.Y() << ", " << nucleon.Z() << ")";

  // make 4-vectors for the kaon and nucleon
  TLorentzVector p4kaon( kaon, sqrt(kaon.Mag2()+mk2) );
  TLorentzVector p4fsnuc( nucleon, sqrt(nucleon.Mag2()+Mf*Mf) );
  // these are in the struck nucleon rest frame...boost them to the lab frame
  p4kaon.Boost( beta );
  p4fsnuc.Boost( beta );

  LOG( "SKHadron", pDEBUG ) <<
    "\nKaon (x,y,z) in lab frame: (" << p4kaon.X() << ", " << p4kaon.Y() << ", " << p4kaon.Z() << ")" <<
    "\nNucleon:                   (" << p4fsnuc.X() << ", " << p4fsnuc.Y() << ", " << p4fsnuc.Z() << ")";

  double pxNf = p4fsnuc.Px();
  double pyNf = p4fsnuc.Py();
  double pzNf = p4fsnuc.Pz();
  double ENf = p4fsnuc.E();

  double pxKf = p4kaon.Px();
  double pyKf = p4kaon.Py();
  double pzKf = p4kaon.Pz();
  double EKf = p4kaon.E();

  //-- Save the particles at the GHEP record

  // mom is mother, not momentum
  int mom = evrec->HitNucleonPosition();

  // AddParticle (int pdg, GHepStatus_t ist, int mom1, int mom2, int dau1, int dau2, double px, double py, double pz, double E, double x, double y, double z, double t)
  
  evrec->AddParticle(
     nuc_pdgc, ist, mom,-1,-1,-1, 
     pxNf, pyNf, pzNf, ENf, vtx.X(), vtx.Y(), vtx.Z(), vtx.T());

  evrec->AddParticle(
     kaon_pdgc,ist, mom,-1,-1,-1, 
     pxKf, pyKf, pzKf, EKf, vtx.X(), vtx.Y(), vtx.Z(), vtx.T());

}

