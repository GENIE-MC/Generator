//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 03, 2007 - CA
   If the formation zone is too large then the particle is placed back,
   just "outside the nucleus". Update the definition of "outside the
   nucleus" to be the maximum distance that a particle can be tracked
   by the intranuclear cascade + a couple of fermis. Brings the code
   in sync with changes in the intranuclear cascade tracking algorithm.
 @ Feb 07, 2009 - CA
   Removed call to AddTargetNucleusRemnant(). This simulation step is now
   performed further upstream in the processing chain.  
 @ Mar 03, 2009 - CA
   Moved into the new DIS package from its previous location (EVGModules).
 @ Sep 15, 2009 - CA
   IsNucleus() is no longer available in GHepParticle. Use pdg::IsIon().
 @ Feb 08, 2013 - CA
   Use the formation zone code from PhysUtils (also used by reweighting) 
   rather than having own implementation here
*/
//____________________________________________________________________________

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Physics/DeepInelastic/EventGen/DISHadronicSystemGenerator.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/Hadronization/HadronizationModelI.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Hadronization/FragmRecUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/PhysUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
using namespace genie::utils;

//___________________________________________________________________________
DISHadronicSystemGenerator::DISHadronicSystemGenerator() :
HadronicSystemGenerator("genie::DISHadronicSystemGenerator")
{

}
//___________________________________________________________________________
DISHadronicSystemGenerator::DISHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::DISHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
DISHadronicSystemGenerator::~DISHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void DISHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  //-- Add an entry for the DIS Pre-Fragm. Hadronic State
  this->AddFinalHadronicSyst(evrec);

  //-- Add the fragmentation products
  this->AddFragmentationProducts(evrec);

  //-- Simulate the formation zone if not taken directly from the 
  //   hadronization model
  this->SimulateFormationZone(evrec);
}
//___________________________________________________________________________
void DISHadronicSystemGenerator::AddFragmentationProducts(
                                                   GHepRecord * evrec) const
{
// Calls a hadronizer and adds the fragmentation products at the GHEP

  GHepParticle * neutrino  = evrec->Probe();
  const TLorentzVector & vtx = *(neutrino->X4());

  //-- Compute the hadronic system invariant mass
  TLorentzVector p4Had = this->Hadronic4pLAB(evrec);
  double W = p4Had.M();

  Interaction * interaction = evrec->Summary();
  interaction->KinePtr()->SetW(W);

  //-- Run the hadronization model and get the fragmentation products:
  //   A collection of ROOT TMCParticles (equiv. to a LUJETS record)

  TClonesArray * plist = fHadronizationModel->Hadronize(interaction);
  if(!plist) {
     LOG("DISHadronicVtx", pWARN) 
                  << "Got an empty particle list. Hadronizer failed!";
     LOG("DISHadronicVtx", pWARN) 
                    << "Quitting the current event generation thread";

     evrec->EventFlags()->SetBitNumber(kHadroSysGenErr, true);

     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Could not simulate the hadronic system");
     exception.SwitchOnFastForward();
     throw exception;

     return;
  }

  //-- Take the hadronic system weight to handle cases that the hadronizer
  //   was asked to produce weighted events
  double wght = fHadronizationModel->Weight();

  //-- Translate the fragmentation products from TMCParticles to
  //   GHepParticles and copy them to the event record.

  int mom = evrec->FinalStateHadronicSystemPosition();
  assert(mom!=-1);
 
  TMCParticle * p = 0;
  TIter particle_iter(plist);

  bool is_nucleus = interaction->InitState().Tgt().IsNucleus();
  GHepStatus_t istfin = (is_nucleus) ?       
                 kIStHadronInTheNucleus : kIStStableFinalState;

  // Vector defining rotation from LAB to LAB' (z:= \vec{phad})
  TVector3 unitvq = p4Had.Vect().Unit();

  // Boost velocity LAB' -> HCM
  TVector3 beta(0,0,p4Had.P()/p4Had.Energy());

  while( (p = (TMCParticle *) particle_iter.Next()) ) {

     int pdgc = p->GetKF();
     int ks   = p->GetKS();

     if(fFilterPreFragmEntries && ks!=1) continue;

     // The fragmentation products are generated in the hadronic CM frame
     // where the z>0 axis is the \vec{phad} direction. For each particle 
     // returned by the hadronizer:
     // - boost it back to LAB' frame {z:=\vec{phad}} / doesn't affect pT
     // - rotate its 3-momentum from LAB' to LAB

     TLorentzVector p4o(p->GetPx(), p->GetPy(), p->GetPz(), p->GetEnergy());
     p4o.Boost(beta); 
     TVector3 p3 = p4o.Vect();
     p3.RotateUz(unitvq); 
     TLorentzVector p4(p3,p4o.Energy());
     
     // copy final state particles to the event record
     GHepStatus_t ist = (ks==1) ? istfin : kIStDISPreFragmHadronicState;

     // handle gammas, and leptons that might come from internal pythia decays
     // mark them as final state particles
     bool not_hadron = (pdgc==kPdgGamma || 
                        pdg::IsNeutralLepton(pdgc) || pdg::IsChargedLepton(pdgc));
     if(not_hadron)  { ist = kIStStableFinalState; }

     int im  = mom + 1 + p->GetParent();
     int ifc = (p->GetFirstChild() == -1) ? -1 : mom + 1 + p->GetFirstChild();
     int ilc = (p->GetLastChild()  == -1) ? -1 : mom + 1 + p->GetLastChild();

     evrec->AddParticle(pdgc, ist, im,-1, ifc, ilc, p4,vtx);

  } // fragmentation-products-iterator

  //-- Handle the case that the hadronizer produced weighted events and
  //   take into account that the current event might be already weighted
  evrec->SetWeight (wght * evrec->Weight());

  plist->Delete();
  delete plist;
}
//___________________________________________________________________________
void DISHadronicSystemGenerator::SimulateFormationZone(
                                                   GHepRecord * evrec) const
{
  LOG("DISHadronicVtx", pDEBUG) 
    << "Simulating formation zone for the DIS hadronic system";
 
  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("DISHadronicVtx", pDEBUG) 
     << "No nuclear target was found - No need to simulate formation zones";
    return;
  }
  // Compute the nuclear radius & how far away a particle is being tracked by
  // the intranuclear hadron transport
  assert(nucltgt && pdg::IsIon(nucltgt->Pdg()));
  double A = nucltgt->A();
  double R = fR0 * TMath::Power(A, 1./3.);
  R *= TMath::Max(fNR,1.); // particle is tracked much further outside the nuclear boundary as the density is non-zero

  // Decay very short living particles so that we can give the formation
  // zone to the daughters
  this->PreHadronTransportDecays(evrec);

  // Get hadronic system's 3-momentum
  GHepParticle * hadronic_system = evrec->FinalStateHadronicSystem();
  TVector3 p3hadr = hadronic_system->P4()->Vect(); // (px,py,pz)

  // Loop over GHEP and set the formation zone to the right particles
  // Limit the maximum formation zone so that particles escaping the
  // nucleus are placed right outside...
  
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  int icurr = -1;
 
  while( (p = (GHepParticle *) piter.Next()) )
  {
    icurr++;

    GHepStatus_t ist = p->Status();
    bool apply_formation_zone = (ist==kIStHadronInTheNucleus);

    if(!apply_formation_zone) continue;

    LOG("DISHadronicVtx", pINFO)
      << "Applying formation-zone to " << p->Name();

    double m = p->Mass();  
    int pdgc = p->Pdg();
    const TLorentzVector & p4 = *(p->P4());
    double ct0=0.;
    pdg::IsNucleon(pdgc) ? ct0=fct0nucleon : ct0=fct0pion; 
    double fz = phys::FormationZone(m,p4,p3hadr,ct0,fK);

    //-- Apply the formation zone step

    double step = fz;
    TVector3 dr = p->P4()->Vect().Unit();            // unit vector along its direction
 // double c  = kLightSpeed / (units::fm/units::ns); // c in fm/nsec
    dr.SetMag(step);                                 // spatial step size
 // double dt = step/c;                              // temporal step:
    double dt = 0;
    TLorentzVector dx4(dr,dt);                       // 4-vector step
    TLorentzVector x4new = *(p->X4()) + dx4;         // new position

    //-- If the formation zone was large enough that the particle is now outside
    //   the nucleus make sure that it is not placed further away from the 
    //   (max distance particles tracked by intranuclear cascade) + ~2 fm
    double epsilon = 2; // fm
    double r       = x4new.Vect().Mag(); // fm
    double rmax    = R+epsilon; 
    if(r > rmax) {
        LOG("DISHadronicVtx", pINFO)
          << "Particle was stepped too far away (r = " << r << " fm)";
        LOG("DISHadronicVtx", pINFO)
          << "Placing it ~2 fm away from the furthermost position tracked "
          << "by intranuclear cascades (r' = " << rmax << " fm)";
        double scale = rmax/r;
        x4new *= scale;
    }

    p->SetPosition(x4new);
  }
}
//___________________________________________________________________________
void DISHadronicSystemGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISHadronicSystemGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISHadronicSystemGenerator::LoadConfig(void)
{
  fHadronizationModel = 0;
  fPreINukeDecayer    = 0;

  //-- Get the requested hadronization model
  fHadronizationModel = 
     dynamic_cast<const HadronizationModelI *> (this->SubAlg("Hadronizer"));
  assert(fHadronizationModel);
 
  //-- Handle pre-intranuke decays
  fPreINukeDecayer =
     dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("PreTransportDecayer"));
  assert(fPreINukeDecayer);

  //-- Get flag to determine whether we copy all fragmentation record entries
  //   into the GHEP record or just the ones marked with kf=1
  GetParamDef( "FilterPreFragm", fFilterPreFragmEntries, false ) ;

  //-- Get parameters controlling the nuclear sizes
  GetParam( "NUCL-R0", fR0 ) ;
  GetParam( "NUCL-NR", fNR ) ;

  //-- Get parameters controlling the formation zone simulation
  GetParam( "FZONE-ct0pion", fct0pion ) ;
  GetParam( "FZONE-ct0nucleon",fct0nucleon ) ;
  GetParam( "FZONE-KPt2", fK ) ;

  LOG("DISHadronicVtx", pDEBUG) << "ct0pion     = " << fct0pion    << " fermi";
  LOG("DISHadronicVtx", pDEBUG) << "ct0nucleon  = " << fct0nucleon << " fermi";
  LOG("DISHadronicVtx", pDEBUG) << "K(pt^2) = " << fK;
}
//____________________________________________________________________________

