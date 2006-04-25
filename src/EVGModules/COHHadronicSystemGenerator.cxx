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

#include <cstdlib>

#include <TVector3.h>

#include "Conventions/Constants.h"
#include "EVGModules/COHHadronicSystemGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/IUtils.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"

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
// This method generates the final state hadronic system (pion + nucleus) in 
// COH interactions
//
  RandomGen * rnd = RandomGen::Instance();

  //-- Access neutrino, initial nucleus and final state prim. lepton entries

  GHepParticle * nu  = evrec->Probe();
  GHepParticle * Ni  = evrec->TargetNucleus();
  GHepParticle * fsl = evrec->FinalStatePrimaryLepton();

  assert(nu);
  assert(Ni);
  assert(fsl);

  //-- Determine the pdg code of the final state pion & nucleus

  int nucl_pdgc = Ni->PdgCode(); // same as the initial nucleus

  Interaction * interaction = evrec->GetInteraction();
  const XclsTag & xcls_tag  = interaction->GetExclusiveTag();

  int pion_pdgc = 0;
  if      (xcls_tag.NPi0()     == 1) pion_pdgc = kPdgPi0;
  else if (xcls_tag.NPiPlus()  == 1) pion_pdgc = kPdgPiPlus;
  else if (xcls_tag.NPiMinus() == 1) pion_pdgc = kPdgPiMinus;
  else {
     LOG("COHHadronicVtx", pFATAL)
                      << "No final state pion information in XclsTag!";
     exit(1);
  }

  //-- Generate the f/s pion and nucleus 4-p

  // get particle masses
  //double Mf   = Ni->Mass(); // Nf:=Ni
  double mpi  = PDGLibrary::Instance()->Find(pion_pdgc)->Mass();
  double mpi2 = TMath::Power(mpi,2);

  // get selected kinematics
  double yo = interaction->GetKinematics().y(true); 
  double to = interaction->GetKinematics().t(true); 

  SLOG("COHHadronicVtx", pDEBUG) << "yo = " << yo << ", to = " << to;

  // 4-momentum transfer q=p(neutrino) - p(f/s lepton)
  TLorentzVector q = *nu->P4() - *fsl->P4();
  TVector3      q3 = q.Vect(); // spatial component

  SLOG("COHHadronicVtx", pDEBUG) 
       << "4-p transfer (q) @ LAB: " << utils::print::P4AsString(&q);

  // compute pion energy and |momentum|
  double Epi  = yo * nu->E();  
  double Epi2 = TMath::Power(Epi,2);
  double ppi2 = Epi2-mpi2;
  double ppi  = TMath::Sqrt(TMath::Max(0.,ppi2));

  SLOG("COHHadronicVtx", pDEBUG)
                      << "f/s pion E = " << Epi << ", |p| = " << ppi;
  assert(Epi>mpi);

  // find angle eta between q and ppi 
  // t = (q4-ppi4)^2 = (q0-ppi0)^2 - (q3-ppi3) 
  // 
  double q0_Epi_2 = TMath::Power(q.E() - Epi, 2);
  double coseta   = (to - q0_Epi_2 + q3.Mag2() + ppi2)/ (2 * q3.Mag() * ppi);
  double sineta   = (1.-TMath::Power(coseta,2));

  SLOG("COHHadronicVtx", pDEBUG) << "cosine(pion, q) = " << coseta;

  // compute transverse and longitudinal ppi components along q
  double ppiL = ppi*coseta;
  double ppiT = ppi*sineta;

  // pion momentum in rotated frame (z along q)
  double phi = 2*kPi* rnd->Random1().Rndm();
  TVector3 ppi3(ppiT*TMath::Cos(phi),ppiT*TMath::Sin(phi),ppiL);

  SLOG("COHHadronicVtx", pDEBUG) 
      << "Pion 3-p with z along q: " << utils::print::Vec3AsString(&ppi3);

  // rotate f/s pion 3-momentum at the Lab

  ppi3.RotateUz(q3.Unit());

  SLOG("COHHadronicVtx", pDEBUG) 
               << "Pion 3-p @ LAB: " << utils::print::Vec3AsString(&ppi3);

  // now figure out the f/s nucleus 4-p

  double pxNf = nu->Px() + Ni->Px() - fsl->Px() - ppi3.Px();
  double pyNf = nu->Py() + Ni->Py() - fsl->Py() - ppi3.Py();
  double pzNf = nu->Pz() + Ni->Pz() - fsl->Pz() - ppi3.Pz();
  double ENf  = nu->E()  + Ni->E()  - fsl->E()  - Epi;

  //-- Save the particles at the GHEP record

  int mom = evrec->TargetNucleusPosition();
  
  evrec->AddParticle(
	    nucl_pdgc,kIStStableFinalState, mom,-1,-1,-1, 
                               	       pxNf, pyNf, pzNf, ENf, 0, 0, 0, 0);
  evrec->AddParticle(
            pion_pdgc,kIStStableFinalState, mom,-1,-1,-1, 
       	                     ppi3.Px(), ppi3.Py(),ppi3.Pz(),Epi, 0,0,0,0);
}
//___________________________________________________________________________

