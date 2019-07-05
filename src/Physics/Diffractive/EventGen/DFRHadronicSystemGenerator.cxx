//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - Feb 15, 2009

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 15, 2009 - CA
   This class was first added in version 2.5.1.

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TVector3.h>

#include "Framework/Conventions/Constants.h"
#include "Physics/Diffractive/EventGen/DFRHadronicSystemGenerator.h"
#include "Framework/EventGen/EVGThreadException.h"
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

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
DFRHadronicSystemGenerator::DFRHadronicSystemGenerator() :
HadronicSystemGenerator("genie::DFRHadronicSystemGenerator")
{

}
//___________________________________________________________________________
DFRHadronicSystemGenerator::DFRHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::DFRHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
DFRHadronicSystemGenerator::~DFRHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void DFRHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system (meson + nucleus) in 
// diffractive scattering 
//
  RandomGen * rnd = RandomGen::Instance();

  //-- Access neutrino, initial nucleus and final state prim. lepton entries
  GHepParticle * nu  = evrec->Probe();
  GHepParticle * Ni  = evrec->HitNucleon();
  GHepParticle * fsl = evrec->FinalStatePrimaryLepton();
  assert(nu);
  assert(Ni);
  assert(fsl);

  Interaction * interaction = evrec->Summary();

  //-- Determine the pdg code of the final state pion 
  const XclsTag & xcls_tag  = interaction->ExclTag();
  int pion_pdgc = 0;
  if      (xcls_tag.NPi0()     == 1) pion_pdgc = kPdgPi0;
  else if (xcls_tag.NPiPlus()  == 1) pion_pdgc = kPdgPiP;
  else if (xcls_tag.NPiMinus() == 1) pion_pdgc = kPdgPiM;
  else {
     LOG("DFRHadronicVtx", pFATAL)
               << "No final state pion information in XclsTag!";
     exit(1);
  }

  //-- Determine the pdg code of the recoil nucleon
  int nucl_pdgc = Ni->Pdg(); // same as the hit  nucleon


  const TLorentzVector & p4nu  = *(nu ->P4());
  const TLorentzVector & p4Ni  = *(Ni ->P4());
  const TLorentzVector & p4fsl = *(fsl->P4());

  // q at LAB
  TLorentzVector q = p4nu - p4fsl;        

  // q at NRF
  TLorentzVector qnrf(q);
  TVector3 beta = p4Ni.BoostVector();
  qnrf.Boost(-beta);

  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  double E    = init_state.ProbeE(kRfHitNucRest);  // neutrino energy
  double M    = target.HitNucMass();
  double mpi  = PDGLibrary::Instance()->Find(pion_pdgc)->Mass();
  double mpi2 = TMath::Power(mpi,2);
  double xo   = interaction->Kine().x(true); 
  double yo   = interaction->Kine().y(true); 
  double to   = interaction->Kine().t(true); 
  double Epi  = qnrf.E() - 0.5*to/M;  
  double Epi2 = TMath::Power(Epi,2);
  double ppi2 = Epi2-mpi2;
  double ppi  = TMath::Sqrt(TMath::Max(0.,ppi2));

  SLOG("DFRHadronicVtx", pINFO) 
         << "Ev = "<< E 
         << ", xo = " << xo << ", yo = " << yo << ", to = " << to;
  SLOG("DFRHadronicVtx", pINFO)
         << "f/s pion E = " << Epi << ", |p| = " << ppi;

  // find angle theta between q and ppi (xi=costheta)
 
  double xi = -to - mpi2 - qnrf.M2() + 2*qnrf.E()*Epi;
  xi /= (2 * TMath::Sqrt(Epi2-mpi2) * qnrf.Vect().Mag());

  if(xi < 0. || xi > 1. || Epi<= mpi) {
     LOG("DFRKinematics", pWARN) << "Invalid selected kinematics; Attempt regenerating";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Invalid selected kinematics");
     exception.SwitchOnStepBack();
     exception.SetReturnStep(0);
     throw exception;
  }

  LOG("DFRHadronicVtx", pINFO) << "|to| = " << -to;
  LOG("DFRHadronicVtx", pINFO) << "q2   = " << qnrf.M2();
  LOG("DFRHadronicVtx", pINFO) << "xi   = " << xi;

  double costheta  = xi;
  double sintheta  = TMath::Sqrt(TMath::Max(0.,1.-xi*xi));

  SLOG("DFRHadronicVtx", pINFO) << "cos(pion, q) = " << costheta;

  // compute transverse and longitudinal ppi components along q
  double ppiL = ppi*costheta;
  double ppiT = ppi*sintheta;
  double phi  = 2*kPi* rnd->RndHadro().Rndm();

  TVector3 ppi3(0,ppiT,ppiL);        // @ NRF' (NRF rotated so that q along z)

  ppi3.RotateZ(phi);                 // randomize transverse components
  ppi3.RotateUz(qnrf.Vect().Unit()); // align longit. component with q in NRF

  TLorentzVector p4pi(ppi3,Epi);
  p4pi.Boost(beta);


  //-- Now figure out the f/s nucleon 4-p

  double pxNf = nu->Px() + Ni->Px() - fsl->Px() - p4pi.Px();
  double pyNf = nu->Py() + Ni->Py() - fsl->Py() - p4pi.Py();
  double pzNf = nu->Pz() + Ni->Pz() - fsl->Pz() - p4pi.Pz();
  double ENf  = nu->E()  + Ni->E()  - fsl->E()  - p4pi.E();

  //-- Save the particles at the GHEP record

  int mom = evrec->HitNucleonPosition();
  const TLorentzVector & vtx = *(nu->X4());

  const Target & tgt = interaction->InitState().Tgt();
  GHepStatus_t ist = (tgt.IsNucleus()) ?
        kIStHadronInTheNucleus : kIStStableFinalState;

  evrec->AddParticle(nucl_pdgc, ist, mom,-1,-1,-1, 
    pxNf, pyNf, pzNf, ENf, 0, 0, 0, 0);

  evrec->AddParticle(pion_pdgc, ist, mom,-1,-1,-1, 
    p4pi.Px(), p4pi.Py(),p4pi.Pz(),p4pi.E(), vtx.X(), vtx.Y(), vtx.Z(), vtx.T());
}
//___________________________________________________________________________

