//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors:

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TVector3.h>

#include "Conventions/Constants.h"
#include "AtharSingleKaon/ASKHadronicSystemGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"
#include "Base/XSecAlgorithmI.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/RunningThreadInfo.h"


using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
ASKHadronicSystemGenerator::ASKHadronicSystemGenerator() :
HadronicSystemGenerator("genie::ASKHadronicSystemGenerator")
{

}
//___________________________________________________________________________
ASKHadronicSystemGenerator::ASKHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::ASKHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
ASKHadronicSystemGenerator::~ASKHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void ASKHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Access cross section algorithm for running thread
  //RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  //const EventGeneratorI * evg = rtinfo->RunningThread();
  //const XSecAlgorithmI *fXSecModel = evg->CrossSectionAlg();
  CalculateHadronicSystem_AtharSingleKaon(evrec);
}
//___________________________________________________________________________
void ASKHadronicSystemGenerator::CalculateHadronicSystem_AtharSingleKaon(GHepRecord * evrec) const
{
//
// This method generates the final state hadronic system (kaon + nucleus) in 
// ASK interactions
//
  //RandomGen * rnd = RandomGen::Instance();

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

  const TLorentzVector & vtx   = *(nu->X4());
  const TLorentzVector & p4nu  = *(nu ->P4());
  const TLorentzVector & p4fsl = *(fsl->P4());

  //-- Determine the pdg code of the final state nucleon
  int nuc_pdgc = (xcls_tag.NProtons()) ? kPdgProton : kPdgNeutron; // there's only ever one nucleon
  int kaon_pdgc = xcls_tag.StrangeHadronPdg();

  const Target & tgt = interaction->InitState().Tgt();
  GHepStatus_t ist = (tgt.IsNucleus()) ?
                          kIStHadronInTheNucleus : kIStStableFinalState;

  //-- basic kinematic inputs
  double M    = (xcls_tag.NProtons()) ? kProtonMass : kNeutronMass; // there's only ever one nucleon
  double mk   = PDGLibrary::Instance()->Find(kaon_pdgc)->Mass(); // K+ and K0 mass are slightly different
  double mk2  = TMath::Power(mk,2);

  //-- specific kinematic quantities
  double kaon_T = kinematics->GetKV(kKVSelTk);
  double kaon_E = kaon_T + mk;
  double pk = sqrt( kaon_E*kaon_E - mk2 );

  TLorentzVector q = p4nu - p4fsl;

  TVector3 qvec = q.Vect(); // this is in lab frame
  double q3 = qvec.Mag(); // in q frame (0,0,q3)

  // Equation 17 of notes from M. Rafi Alam dated 6 November 2013
  double eN = q.E() + M - kaon_E; // nucleon total energy
  double cos_thetaKq = (q3*q3 + pk*pk + M*M - eN*eN)/(2*q3*pk);
  // this can be slightly larger than 1 due to numerical precision issues -- don't let it be
  if( cos_thetaKq > 1.0 ) cos_thetaKq = 1.0;

  // Get phi for the k-q plane relative to nu-l plane
  double phi_kq = kinematics->GetKV(kKVSelphikq);

  TVector3 kaon;
  kaon.SetMagThetaPhi( pk, TMath::ACos(cos_thetaKq), phi_kq );

  TVector3 nucleon( -kaon.X(), -kaon.Y(), q3 - kaon.Z() ); // force 3-momentum conservation in q frame

  LOG( "ASKHadron", pDEBUG ) <<
    "\nTk " << kaon_T << " Tlep " << kinematics->GetKV(kKVSelTl) << " costhetalep " << kinematics->GetKV(kKVSelctl) << " philep " << p4fsl.Vect().Phi() << " phikq " << phi_kq <<
    "\nKaon (x,y,z) in q frame: (" << kaon.X() << ", " << kaon.Y() << ", " << kaon.Z() << ")" <<
    "\nNucleon:                 (" << nucleon.X() << ", " << nucleon.Y() << ", " << nucleon.Z() << ")";

  // Transform hadron momenta from z axis along q to lab frame
  kaon.RotateUz(qvec.Unit());
  nucleon.RotateUz(qvec.Unit());

  LOG( "ASKHadron", pDEBUG ) <<
    "\nKaon (x,y,z) in lab frame: (" << kaon.X() << ", " << kaon.Y() << ", " << kaon.Z() << ")" <<
    "\nNucleon:                   (" << nucleon.X() << ", " << nucleon.Y() << ", " << nucleon.Z() << ")";

  double pxNf = nucleon.Px();
  double pyNf = nucleon.Py();
  double pzNf = nucleon.Pz();
  double ENf = sqrt(pxNf*pxNf + pyNf*pyNf + pzNf*pzNf + M*M);  

  double pxKf = kaon.Px();
  double pyKf = kaon.Py();
  double pzKf = kaon.Pz();
  double EKf = sqrt(pxKf*pxKf + pyKf*pyKf + pzKf*pzKf + mk2); 

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

