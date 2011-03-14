//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Sep 22, 2008

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 22, 2008 - CA
   This event generation modules was first added in version 2.5.1 as part of
   the new event generation thread handling MEC interactions. 

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/RunningThreadInfo.h"
#include "EVGCore/EventGeneratorI.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepFlags.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "MEC/MECGenerator.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;

//___________________________________________________________________________
MECGenerator::MECGenerator() :
EventRecordVisitorI("genie::MECGenerator")
{

}
//___________________________________________________________________________
MECGenerator::MECGenerator(string config) :
EventRecordVisitorI("genie::MECGenerator", config)
{

}
//___________________________________________________________________________
MECGenerator::~MECGenerator()
{

}
//___________________________________________________________________________
void MECGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  this -> SelectKinematics    (evrec);
  this -> AddNucleonCluster   (evrec);
  this -> AddTargetRemnant    (evrec);
  this -> AddFinalStateLepton (evrec);
  this -> DecayNucleonCluster (evrec);

  LOG("MEC", pNOTICE) << *evrec;
}
//___________________________________________________________________________
void MECGenerator::SelectKinematics(GHepRecord * evrec) const
{
  // Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  Interaction * interaction = evrec->Summary();

  // Random num generator
  RandomGen * rnd = RandomGen::Instance();

  // Hardcode bogus limits for the time-being
  double Q2min =  0.01;
  double Q2max = 10.00;
  double Wmin  =  0.50;
  double Wmax  =  1.50;

  // Scan for maximum differential cross section at current energy
  const int nq=30;
  const int nw=20;
  double dQ2 = (Q2max-Q2min) / (nq-1);
  double dW  = (Wmax-Wmin )  / (nw-1);
  double xsec_max =  0;
  for(int iw=0; iw<nw; iw++) {
    for(int iq=0; iq<nq; iq++) {
      double Q2 = Q2min + iq*dQ2;
      double W  = Wmin  + iw*dW;
      interaction->KinePtr()->SetQ2(Q2);  
      interaction->KinePtr()->SetW (W);   
      double xsec = fXSecModel->XSec(interaction, kPSWQ2fE);
      xsec_max = TMath::Max(xsec, xsec_max);
    }
  }
  LOG("MEC", pNOTICE) << "xsec_max = " << xsec_max;

  // Select kinematics
  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("MEC", pWARN)
           << "Couldn't select a valid W, Q^2 pair after " 
           << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     // Generate next pair
     double gQ2 = Q2min + (Q2max-Q2min) * rnd->RndKine().Rndm();
     double gW  = Wmin  + (Wmax -Wmin ) * rnd->RndKine().Rndm();

     // Calculate d2sigma/dQ2dW
     interaction->KinePtr()->SetQ2(gQ2);  
     interaction->KinePtr()->SetW (gW);   
     double xsec = fXSecModel->XSec(interaction, kPSWQ2fE);
     
     // Decide whether to accept the current kinematics
     double t = xsec_max * rnd->RndKine().Rndm();
     double J = 1; // jacobean
     accept = (t < J*xsec);

     // If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("MEC", pINFO) << "Selected: Q^2 = " << gQ2 << ", W = " << gW;

        double gx = 0;
        double gy = 0;
        double E  = evrec->Probe()->E();
        kinematics::WQ2toXY(E,kNucleonMass,gW,gQ2,gx,gy);

        // lock selected kinematics & clear running values
        interaction->KinePtr()->SetQ2(gQ2, true);
        interaction->KinePtr()->SetW (gW,  true);
        interaction->KinePtr()->Setx (gx,  true);
        interaction->KinePtr()->Sety (gy,  true);
        interaction->KinePtr()->ClearRunningValues();
        
        return;
     }//accepted?
  }//iter
}
//___________________________________________________________________________
void MECGenerator::AddFinalStateLepton(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  // Look-up selected kinematics
  double Q2 = interaction->Kine().Q2(true);
  double y  = interaction->Kine().y(true);

  // Auxiliary params
  double Ev  = init_state.ProbeE(kRfLab);
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  double El  = (1-y)*Ev;
  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)

  LOG("MEC", pNOTICE)
          << "fsl: E = " << El << ", |p//| = " << plp << "[pT] = " << plt;

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);

  // Take a unit vector along the neutrino direction
  TVector3 unit_nudir = evrec->Probe()->P4()->Vect().Unit();

  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the LAB
  TVector3 p3l(pltx,plty,plp);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in LAB
  TLorentzVector p4l(p3l,El);

  // Figure out the Final State Lepton PDG Code
  int pdgc = interaction->FSPrimLepton()->PdgCode();

  // Lepton 4-position (= interacton vtx)
  TLorentzVector v4(*evrec->Probe()->X4());

  evrec->AddParticle(pdgc, kIStStableFinalState, -1, 0, -1, -1, p4l, v4);
}
//___________________________________________________________________________
void MECGenerator::AddNucleonCluster(GHepRecord * evrec) const
{
  const int nc = 3;

  int    clusters [nc] = { kPdgClusterNN, kPdgClusterNP, kPdgClusterPP };
  double prob     [nc] = { 0.25,          0.50,          0.25          };

  GHepParticle * target = evrec->TargetNucleus();
  int Z = target->Z();
  int A = target->A();
  int np = Z;
  int nn = A-Z;

  RandomGen * rnd = RandomGen::Instance();

  bool selected = false;
  int  cluster  = 0;

  while(!selected) {
     double prob_gen = rnd->RndGen().Rndm();
     double prob_sum = 0;
     for(int ic = 0; ic < nc; ic++) {
        prob_sum += prob[ic];
        if(prob_gen < prob_sum) {
           cluster = clusters[ic];
           break;
        }
     }
     if ( cluster==kPdgClusterNN && nn>=2          ) selected=true;
     if ( cluster==kPdgClusterNP && nn>=1 && np>=1 ) selected=true;
     if ( cluster==kPdgClusterPP && np>=2          ) selected=true;
  }

  double mc = PDGLibrary::Instance()->Find(cluster)->Mass();

  const TLorentzVector p4(0.,0.,0., mc);
  const TLorentzVector v4(0.,0.,0., 0.);

  LOG("MEC", pINFO) << "Adding nucleon cluster [pdgc = " << cluster << "]";

  evrec->AddParticle(cluster, kIStNucleonTarget, -1, -1, -1, -1, p4, v4);
}
//___________________________________________________________________________
void MECGenerator::AddTargetRemnant(GHepRecord * evrec) const
{
  GHepParticle * target  = evrec->TargetNucleus();
  GHepParticle * cluster = evrec->Particle(2);

  int Z = target->Z();
  int A = target->A();

  if(cluster->Pdg() == kPdgClusterNN) { A-=2; ;     }
  if(cluster->Pdg() == kPdgClusterNP) { A-=2; Z-=1; }
  if(cluster->Pdg() == kPdgClusterPP) { A-=2; Z-=2; }

  int ipdgc = pdg::IonPdgCode(A, Z);

  const TLorentzVector & p4cluster = *(cluster->P4());
  const TLorentzVector & p4tgt     = *(target ->P4());

  const TLorentzVector p4 = p4tgt - p4cluster;
  const TLorentzVector v4(0.,0.,0., 0.);

  int mom = evrec->TargetNucleusPosition();
  evrec->AddParticle(ipdgc,kIStStableFinalState, mom,-1,-1,-1, p4,v4);  
}
//___________________________________________________________________________
void MECGenerator::DecayNucleonCluster(GHepRecord * evrec) const
{
  int cluster_pos = 2;

  GHepParticle * cluster = evrec->Particle(cluster_pos);

  int nuc1=0, nuc2=0;

  if(cluster->Pdg() == kPdgClusterNN) { nuc1=kPdgNeutron; nuc2=kPdgNeutron; }
  if(cluster->Pdg() == kPdgClusterNP) { nuc1=kPdgNeutron; nuc2=kPdgProton;  }
  if(cluster->Pdg() == kPdgClusterPP) { nuc1=kPdgProton;  nuc2=kPdgProton;  }

  int mom = cluster_pos;

  double m1 = PDGLibrary::Instance()->Find(nuc1)->Mass();
  double m2 = PDGLibrary::Instance()->Find(nuc2)->Mass();

  const TLorentzVector p4_1(0.,0., 1., TMath::Sqrt(1.+m1*m1));
  const TLorentzVector p4_2(0.,0.,-1., TMath::Sqrt(1.+m2*m2));

  const TLorentzVector v4(0.,0.,0.,0.);

  evrec->AddParticle(nuc1, kIStHadronInTheNucleus, mom,-1,-1,-1, p4_1, v4);  
  evrec->AddParticle(nuc2, kIStHadronInTheNucleus, mom,-1,-1,-1, p4_2, v4);  
}
//___________________________________________________________________________
void MECGenerator::Configure(const Registry & config)   
{
  Algorithm::Configure(config);
  this->LoadConfig();
} 
//___________________________________________________________________________ 
void MECGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void MECGenerator::LoadConfig(void)
{

}
//___________________________________________________________________________

