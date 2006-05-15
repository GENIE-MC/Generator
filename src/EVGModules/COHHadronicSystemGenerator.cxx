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

  Interaction * interaction = evrec->GetInteraction();
  const XclsTag & xcls_tag  = interaction->GetExclusiveTag();

  //-- Access neutrino, initial nucleus and final state prim. lepton entries
  GHepParticle * nu  = evrec->Probe();
  GHepParticle * Ni  = evrec->TargetNucleus();
  GHepParticle * fsl = evrec->FinalStatePrimaryLepton();
  assert(nu);
  assert(Ni);
  assert(fsl);

  const TLorentzVector & p4nu  = *(nu ->P4());
  const TLorentzVector & p4fsl = *(fsl->P4());

  //-- Determine the pdg code of the final state pion & nucleus
  int nucl_pdgc = Ni->PdgCode(); // same as the initial nucleus
  int pion_pdgc = 0;
  if      (xcls_tag.NPi0()     == 1) pion_pdgc = kPdgPi0;
  else if (xcls_tag.NPiPlus()  == 1) pion_pdgc = kPdgPiPlus;
  else if (xcls_tag.NPiMinus() == 1) pion_pdgc = kPdgPiMinus;
  else {
     LOG("COHHadronicVtx", pFATAL)
               << "No final state pion information in XclsTag!";
     exit(1);
  }

  //-- basic kinematic inputs
  double E    = nu->E();  
  double M    = kNucleonMass;
  double mpi  = PDGLibrary::Instance()->Find(pion_pdgc)->Mass();
  double mpi2 = TMath::Power(mpi,2);
  double xo   = interaction->GetKinematics().x(true); 
  double yo   = interaction->GetKinematics().y(true); 
  double to   = interaction->GetKinematics().t(true); 

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
  //   the outgoing pion  Rein & Seghal, Nucl.Phys.B223.29-44(1983), p.35

  TLorentzVector q = p4nu - p4fsl;

  SLOG("COHHadronicVtx", pINFO) 
          << "\n 4-p transfer q @ LAB: " << utils::print::P4AsString(&q);

  //-- find angle theta between q and ppi (xi=costheta)
  //   note: t=|(ppi-q)^2|, Rein & Seghal, Nucl.Phys.B223.29-44(1983), p.36
 
  double xi = 1. + M*xo/Epi - 0.5*mpi2/Epi2 - 0.5*to/Epi2;
  xi /= TMath::Sqrt((1.+2.*M*xo/Epi)*(1.-mpi2/Epi2));

  double costheta  = xi;
  double sintheta  = TMath::Sqrt(TMath::Max(0.,1.-xi*xi));

  SLOG("COHHadronicVtx", pINFO) << "cos(pion, q) = " << costheta;

  // compute transverse and longitudinal ppi components along q
  double ppiL = ppi*costheta;
  double ppiT = ppi*sintheta;

  double phi = 2*kPi* rnd->Random1().Rndm();

  TVector3 ppi3(0,ppiT,ppiL);

  ppi3.RotateUz(q.Vect().Unit()); // align longit. component with q in LAB
  ppi3.RotateZ(phi);              // randomize transverse components

  SLOG("COHHadronicVtx", pINFO) 
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

