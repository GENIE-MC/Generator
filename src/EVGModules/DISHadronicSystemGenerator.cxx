//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/DISHadronicSystemGenerator.h"
#include "Fragmentation/HadronizationModelI.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/FragmRecUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
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

  //-- If the struck nucleon was within a nucleus, then add the final state
  //   nucleus at the EventRecord
  this->AddTargetNucleusRemnant(evrec);

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
 
  // Get hadronic system's 3-momentum
  GHepParticle * hadronic_system = evrec->FinalStateHadronicSystem();
  TVector3 p3hadr = hadronic_system->P4()->Vect(); // (px,py,pz)

  // Loop over GHEP and run intranuclear rescattering on handled particles
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  int icurr = -1;
 
  while( (p = (GHepParticle *) piter.Next()) )
  {
    icurr++;

   //-- Decide whether we need to apply the formation zone
   //-- Applied to direct descendants of the the 'HadronicSystem' GHEP entry 
   //-- -in case od the KNO model- or descedants of  the JETSET special particles 
   //-- -cluster,string,indep-)
      
   bool apply_formation_zone = false;
 
   // mom 
   int imom = p->FirstMother();
   if(imom<0) continue;     

   // grand-mom pdgc
   int mom_pdg = evrec->Particle(imom)->Pdg();
    
   if (mom_pdg == kPdgHadronicSyst ||
       mom_pdg == kPdgCluster      ||
       mom_pdg == kPdgString       ||
       mom_pdg == kPdgIndep) apply_formation_zone = true;
   if(!apply_formation_zone) continue;

    //-- Compute the formation zone
    //
    TVector3 p3  = p->P4()->Vect();      // hadron's: p (px,py,pz)
    double   m   = p->Mass();            //           m
    double   m2  = m*m;                  //           m^2
    double   P   = p->P4()->P();         //           |p|
    double   Pt  = p3.Pt(p3hadr);        //           pT
    double   Pt2 = Pt*Pt;                //           pT^2
    double   fz  = P*fct0*m/(m2+fK*Pt2); //           formation zone, in m

    LOG("DISHadronicVtx", pNOTICE)
      << p->Name() << ": |P| = " << P << " GeV, Pt = " << Pt
                                << " GeV, Formation Zone = " << fz << " m";

    //-- Apply the formation zone

    double step = fz;

    TVector3 dr = p->P4()->Vect().Unit();          // unit vector along its direction
    double c  = kLightSpeed / (units::m/units::s); // c in m/sec
    dr.SetMag(step);                               // spatial step size
    double dt = step/c;                            // temporal step:
    TLorentzVector dx4(dr,dt);                     // 4-vector step
    TLorentzVector x4new = *(p->X4()) + dx4;       // new position

    LOG("DISHadronicVtx", pDEBUG)
         << "\n Init direction = " << print::Vec3AsString(&dr)
         << "\n Init position (in m,sec) = " << print::X4AsString(p->X4())
         << "\n Fin  position (in m,sec) = " << print::X4AsString(&x4new);

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
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fHadronizationModel = 0;

  //-- Get the requested hadronization model
  fHadronizationModel = 
     dynamic_cast<const HadronizationModelI *> (this->SubAlg("Hadronizer"));

  assert(fHadronizationModel);

  //-- Get flag to determine whether we copy all fragmentation record entries
  //   into the GHEP record or just the ones marked with kf=1
  fFilterPreFragmEntries = fConfig->GetBoolDef("FilterPreFragm",false);

  //-- Get parameters controlling the formation zone simulation
  //
  fct0 = fConfig->GetDoubleDef ("ct0",  gc->GetDouble("FZONE-ct0")); // fm
  fK   = fConfig->GetDoubleDef ("Kpt2", gc->GetDouble("FZONE-KPt2"));

  LOG("DISHadronicVtx", pDEBUG) << "ct0     = " << fct0 << " fermi";
  LOG("DISHadronicVtx", pDEBUG) << "K(pt^2) = " << fK;
}
//____________________________________________________________________________

