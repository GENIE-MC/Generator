//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the new Coherent package from its previous location (EVGModules 
   package)
 @ Mar 03, 2009 - CA
   Renamed COHPiHadronicSystemGenerator -> COHHadronicSystemGenerator in
   anticipation of reusing the code for simulating coherent production of
   vector mesons.
 @ Apr 02, 2009 - CA,HG,PK
   Bug fix: Reverse the order of the pion momentum rotations: Randomize the
   transverse component direction in the x'y' plane before aligning z' with 
   the direction of the momentum transfer q in the LAB.
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TVector3.h>

#include "Framework/Conventions/Constants.h"
#include "Physics/Coherent/EventGen/COHHadronicSystemGenerator.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
COHHadronicSystemGenerator::COHHadronicSystemGenerator() :
  HadronicSystemGenerator("genie::COHHadronicSystemGenerator")
{

}
//___________________________________________________________________________
COHHadronicSystemGenerator::COHHadronicSystemGenerator(string config) :
  HadronicSystemGenerator("genie::COHHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
COHHadronicSystemGenerator::~COHHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void COHHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  // Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  const XSecAlgorithmI *fXSecModel = evg->CrossSectionAlg();
  if (fXSecModel->Id().Name() == "genie::ReinSehgalCOHPiPXSec") {
    CalculateHadronicSystem_ReinSehgal(evrec);
  } else if ((fXSecModel->Id().Name() == "genie::BergerSehgalCOHPiPXSec2015")) {
    CalculateHadronicSystem_BergerSehgal(evrec);
  } else if ((fXSecModel->Id().Name() == "genie::BergerSehgalFMCOHPiPXSec2015")) {
    CalculateHadronicSystem_BergerSehgalFM(evrec);
  } else if ((fXSecModel->Id().Name() == "genie::AlvarezRusoCOHPiPXSec")) {
    CalculateHadronicSystem_AlvarezRuso(evrec);
  }
  else {
    LOG("COHHadronicSystemGenerator",pFATAL) <<
      "ProcessEventRecord >> Cannot calculate hadronic system for " <<
      fXSecModel->Id().Name();
  }
}
//___________________________________________________________________________
int COHHadronicSystemGenerator::getPionPDGCodeFromXclTag(const XclsTag& xcls_tag) const
{
  int pion_pdgc = 0;
  if      (xcls_tag.NPi0()     == 1) pion_pdgc = kPdgPi0;
  else if (xcls_tag.NPiPlus()  == 1) pion_pdgc = kPdgPiP;
  else if (xcls_tag.NPiMinus() == 1) pion_pdgc = kPdgPiM;
  else {
    LOG("COHHadronicVtx", pFATAL)
      << "No final state pion information in XclsTag!";
    exit(1);
  }
  return pion_pdgc;
}
//___________________________________________________________________________
void COHHadronicSystemGenerator::CalculateHadronicSystem_BergerSehgal(GHepRecord * evrec) const
{
  // Treatment of the hadronic side is identical to Rein-Sehgal if we assume an infinite 
  // mass for the nucleus.
  CalculateHadronicSystem_ReinSehgal(evrec);
}
//___________________________________________________________________________
void COHHadronicSystemGenerator::CalculateHadronicSystem_BergerSehgalFM(GHepRecord * evrec) const
{
  //
  // This method generates the final state hadronic system (pion + nucleus) in 
  // COH interactions
  //
  RandomGen * rnd = RandomGen::Instance();

  Interaction * interaction = evrec->Summary();
  const XclsTag & xcls_tag  = interaction->ExclTag();
  const InitialState & init_state = interaction -> InitState();

  //-- Access neutrino, initial nucleus and final state prim. lepton entries
  GHepParticle * nu  = evrec->Probe();
  GHepParticle * Ni  = evrec->TargetNucleus();
  GHepParticle * fsl = evrec->FinalStatePrimaryLepton();
  assert(nu);
  assert(Ni);
  assert(fsl);

  const TLorentzVector & vtx   = *(nu->X4());
  const TLorentzVector & p4nu  = *(nu ->P4());
  const TLorentzVector & p4fsl = *(fsl->P4());

  //-- Determine the pdg code of the final state pion & nucleus
  int nucl_pdgc = Ni->Pdg(); // same as the initial nucleus
  int pion_pdgc = getPionPDGCodeFromXclTag(xcls_tag);

  //-- basic kinematic inputs
  double E    = nu->E();  
  double Q2   = interaction->Kine().Q2(true);
  double y    = interaction->Kine().y(true); 
  double t    = interaction->Kine().t(true); 
  double MA   = init_state.Tgt().Mass(); 
  // double MA2  = TMath::Power(MA, 2.);   // Unused
  double mpi  = PDGLibrary::Instance()->Find(pion_pdgc)->Mass();
  double mpi2 = TMath::Power(mpi,2);

  SLOG("COHHadronicVtx", pINFO) 
    << "Ev = "<< E << ", Q^2 = " << Q2 
    << ", y = " << y << ", t = " << t;

  double Epi = y * E - t / (2 * MA);
  double ppi2   = Epi * Epi - mpi2;
  double ppi    = ppi2 > 0.0 ? TMath::Sqrt(ppi2) : 0.0;

  double costheta = (t - Q2 - mpi2) / (2 * ( (y *E - Epi) * Epi - 
       ppi * sqrt(TMath::Power(y * E - Epi, 2.) + t) ) );

  if ((costheta > 1.0) || (costheta < -1.0)) {
    SLOG("COHHadronicVtx", pERROR) 
      << "Unphysical pion angle!";
  }

  double sintheta = TMath::Sqrt(1 - costheta * costheta);

  //-- first work in the c.m.s. frame
  // double S        = 2 * MA * nuh - Q2 + MA2;
  // double S_2      = S >= 0 ? TMath::Sqrt(S) : 0.0;  // TODO - Error here?
  // double Pcm      = MA * TMath::Sqrt( (nuh*nuh + Q2)/S );
  // double Epi      = (S + mpi2 - MA2)/(2 * S_2);
  // double EAprime  = (S - mpi2 + MA2)/(2 * S_2);
  // double EA       = (S + MA2 + Q2)/(2 * S_2);
  // double PAprime2 = TMath::Power(EAprime,2.0) - MA2;
  // double PAprime  = TMath::Sqrt(PAprime2);
  // double tA       = TMath::Power((EAprime - EA),2.0) - TMath::Power(PAprime,2.0) - 
  //   TMath::Power(Pcm, 2.0);
  // double tB       = 2 * Pcm * PAprime;
  // double cosT     = (t - tA)/tB;
  // double sinT     = TMath::Sqrt(1 - cosT*cosT);
  // double PAz      = PAprime * cosT;
  // double PAperp   = PAprime * sinT;
  // double PPiz     = -PAz;

  // Randomize transverse components
  double phi    = 2 * kPi * rnd->RndHadro().Rndm();
  double ppix   = ppi * sintheta * TMath::Cos(phi);
  double ppiy   = ppi * sintheta * TMath::Sin(phi);
  double ppiz   = ppi * costheta;

  // boost back to the lab frame
  // double beta      = TMath::Sqrt( nuh*nuh + Q2 )/(nuh + MA);
  // double gamma     = (nuh + MA)/TMath::Sqrt(S);
  // double betagamma = beta * gamma;

  // double epi  = gamma*Epi + betagamma*PPiz;
  // double ppiz = betagamma*Epi + gamma*PPiz;

  // double ea  = gamma*EAprime + betagamma*PAz;
  // double paz = betagamma*EAprime + gamma*PAz;

  // Now rotate so our axes are aligned with the lab instead of q
  TLorentzVector q = p4nu - p4fsl;
  TVector3 ppi3(ppix, ppiy, ppiz);
  ppi3.RotateUz(q.Vect().Unit());

  // Nucleus...
  // TVector3 pa(PAx,PAy,paz);
  // pa.RotateUz(q.Vect().Unit());

  // now figure out the f/s nucleus 4-p

  double pxNf = nu->Px() + Ni->Px() - fsl->Px() - ppi3.Px();
  double pyNf = nu->Py() + Ni->Py() - fsl->Py() - ppi3.Py();
  double pzNf = nu->Pz() + Ni->Pz() - fsl->Pz() - ppi3.Pz();
  double ENf  = nu->E()  + Ni->E()  - fsl->E()  - Epi;

  //-- Save the particles at the GHEP record

  int mom = evrec->TargetNucleusPosition();

  // Nucleus - need to balance overall 4-momentum
  evrec->AddParticle(nucl_pdgc,kIStStableFinalState, mom,-1,-1,-1, 
                     pxNf, pyNf, pzNf, ENf, 0, 0, 0, 0);

  // evrec->AddParticle(
  //     nucl_pdgc,kIStStableFinalState, mom,-1,-1,-1,
  //     pa.Px(), pa.Py(), pa.Pz(), ea, 0, 0, 0, 0);

  evrec->AddParticle(
      pion_pdgc,kIStStableFinalState, mom,-1,-1,-1,
      ppi3.Px(), ppi3.Py(), ppi3.Pz(), Epi, vtx.X(), vtx.Y(), vtx.Z(), vtx.T());
}
//___________________________________________________________________________
void COHHadronicSystemGenerator::CalculateHadronicSystem_ReinSehgal(GHepRecord * evrec) const
{
  //
  // This method generates the final state hadronic system (pion + nucleus) in 
  // COH interactions
  //
  RandomGen * rnd = RandomGen::Instance();

  Interaction * interaction = evrec->Summary();
  const XclsTag & xcls_tag  = interaction->ExclTag();

  //-- Access neutrino, initial nucleus and final state prim. lepton entries
  GHepParticle * nu  = evrec->Probe();
  GHepParticle * Ni  = evrec->TargetNucleus();
  GHepParticle * fsl = evrec->FinalStatePrimaryLepton();
  assert(nu);
  assert(Ni);
  assert(fsl);

  const TLorentzVector & vtx   = *(nu->X4());
  const TLorentzVector & p4nu  = *(nu ->P4());
  const TLorentzVector & p4fsl = *(fsl->P4());

  //-- Determine the pdg code of the final state pion & nucleus
  int nucl_pdgc = Ni->Pdg(); // same as the initial nucleus
  int pion_pdgc = getPionPDGCodeFromXclTag(xcls_tag);

  //-- basic kinematic inputs
  double E    = nu->E();  
  double M    = kNucleonMass;
  double mpi  = PDGLibrary::Instance()->Find(pion_pdgc)->Mass();
  double mpi2 = TMath::Power(mpi,2);
  double xo   = interaction->Kine().x(true); 
  double yo   = interaction->Kine().y(true); 
  double to   = interaction->Kine().t(true); 

  SLOG("COHHadronicVtx", pINFO) 
    << "Ev = "<< E << ", xo = " << xo 
    << ", yo = " << yo << ", to = " << to;

  //-- compute pion energy and |momentum|
  double Epi  = yo * E;  
  double Epi2 = TMath::Power(Epi,2);
  double ppi2 = Epi2-mpi2;
  double ppi  = TMath::Sqrt(TMath::Max(0.,ppi2));

  SLOG("COHHadronicVtx", pINFO)
    << "f/s pion E = " << Epi << ", |p| = " << ppi;
  assert(Epi>mpi);

  //-- 4-momentum transfer q=p(neutrino) - p(f/s lepton)
  //   Note: m^2 = q^2 < 0
  //   Also, since the nucleus is heavy all energy loss is comminicated to
  //   the outgoing pion  Rein & Sehgal, Nucl.Phys.B223.29-44(1983), p.35

  TLorentzVector q = p4nu - p4fsl;

  SLOG("COHHadronicVtx", pINFO) 
    << "\n 4-p transfer q @ LAB: " << utils::print::P4AsString(&q);

  //-- find angle theta between q and ppi (xi=costheta)
  //   note: t=|(ppi-q)^2|, Rein & Sehgal, Nucl.Phys.B223.29-44(1983), p.36

  double xi = 1. + M*xo/Epi - 0.5*mpi2/Epi2 - 0.5*to/Epi2;
  xi /= TMath::Sqrt((1.+2.*M*xo/Epi)*(1.-mpi2/Epi2));

  double costheta  = xi;
  double sintheta  = TMath::Sqrt(TMath::Max(0.,1.-xi*xi));

  SLOG("COHHadronicVtx", pINFO) << "cos(pion, q) = " << costheta;

  // compute transverse and longitudinal ppi components along q
  double ppiL = ppi*costheta;
  double ppiT = ppi*sintheta;

  double phi = 2*kPi* rnd->RndHadro().Rndm();

  TVector3 ppi3(0,ppiT,ppiL);

  ppi3.RotateZ(phi);              // randomize transverse components
  ppi3.RotateUz(q.Vect().Unit()); // align longit. component with q in LAB

  SLOG("COHHadronicVtx", pINFO) 
    << "Pion 3-p @ LAB: " << utils::print::Vec3AsString(&ppi3);

  // now figure out the f/s nucleus 4-p

  double pxNf = nu->Px() + Ni->Px() - fsl->Px() - ppi3.Px();
  double pyNf = nu->Py() + Ni->Py() - fsl->Py() - ppi3.Py();
  double pzNf = nu->Pz() + Ni->Pz() - fsl->Pz() - ppi3.Pz();
  double ENf  = nu->E()  + Ni->E()  - fsl->E()  - Epi;

  //-- Save the particles at the GHEP record

  int mom = evrec->TargetNucleusPosition();

  evrec->AddParticle(nucl_pdgc,kIStStableFinalState, mom,-1,-1,-1, 
                     pxNf, pyNf, pzNf, ENf, 0, 0, 0, 0);

  evrec->AddParticle(pion_pdgc,kIStStableFinalState, mom,-1,-1,-1, 
                     ppi3.Px(), ppi3.Py(),ppi3.Pz(),Epi, vtx.X(), vtx.Y(), vtx.Z(), vtx.T());
}
//___________________________________________________________________________
void COHHadronicSystemGenerator::CalculateHadronicSystem_AlvarezRuso(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->Summary();
  const Kinematics &   kinematics = interaction -> Kine();
  GHepParticle * nu  = evrec->Probe();
  GHepParticle * Ni  = evrec->TargetNucleus();
  GHepParticle * fsl = evrec->FinalStatePrimaryLepton();

  // Pion
  const TLorentzVector ppi  = kinematics.HadSystP4();
  const TVector3 ppi3 = ppi.Vect();
  const double Epi = ppi.E();
  int pion_pdgc=0;
  if ( interaction->ProcInfo().IsWeakCC() ) {
    if( nu->Pdg() > 0 ){ // neutrino
      pion_pdgc = kPdgPiP;
    }
    else{ // anti-neutrino
      pion_pdgc = kPdgPiM;
    }
  }
  else if ( interaction->ProcInfo().IsWeakNC() ) {
    pion_pdgc = kPdgPi0;
  }
  else{
    LOG("COHHadronicSystemGeneratorAR", pFATAL)
      << "Could not determine pion involved in interaction";
    exit(1);
  }

  //
  // Nucleus
  int nucl_pdgc = Ni->Pdg(); // pdg of final nucleus same as the initial nucleus
  double pxNf = nu->Px() + Ni->Px() - fsl->Px() - ppi3.Px();
  double pyNf = nu->Py() + Ni->Py() - fsl->Py() - ppi3.Py();
  double pzNf = nu->Pz() + Ni->Pz() - fsl->Pz() - ppi3.Pz();
  double ENf  = nu->E()  + Ni->E()  - fsl->E()  - Epi;
  //
  // Both
  const TLorentzVector & vtx   = *(nu->X4());
  int mom = evrec->TargetNucleusPosition();

  //
  // Fill the records
  evrec->AddParticle(nucl_pdgc,kIStStableFinalState, mom,-1,-1,-1,
                     pxNf, pyNf, pzNf, ENf, 0, 0, 0, 0);

  evrec->AddParticle(pion_pdgc,kIStStableFinalState, mom,-1,-1,-1,
                     ppi3.Px(), ppi3.Py(),ppi3.Pz(),Epi, vtx.X(), vtx.Y(), vtx.Z(), vtx.T());
}

