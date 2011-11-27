//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory               

         Steve Dytman <dytman+ \at pitt.edu>
         Pittsburgh University

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 22, 2008 - CA
   Skeleton was first added in version 2.5.1
 @ Nov 24, 2010 - CA

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
#include "Nuclear/NuclearModelI.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"
#include "Utils/PrintUtils.h"

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
void MECGenerator::ProcessEventRecord(GHepRecord * event) const
{
  this -> AddNucleonCluster     (event);
  this -> AddTargetRemnant      (event);
  this -> GenerateFermiMomentum (event);
  this -> SelectKinematics      (event);
  this -> AddFinalStateLepton   (event);
  this -> RecoilNucleonCluster  (event);
  this -> DecayNucleonCluster   (event);

  LOG("MEC", pNOTICE) << *event;
}
//___________________________________________________________________________
void MECGenerator::AddNucleonCluster(GHepRecord * event) const
{
// Add a p+p, p+n or n+n nucleon cluster as the "hit" object.
// Using built-in probabilities for each nucleon cluster.
// The function checks that there is a sufficient number of neutrons and
// protons in the hit nucleus for the selected nucleon cluster.
//
  const int nc = 3;
  int    clusters [nc] = { kPdgClusterNN, kPdgClusterNP, kPdgClusterPP };
  double prob     [nc] = { 0.12,          0.88,          0.12          };
  int nupdg = event->Probe()->Pdg();
  if (pdg::IsNeutrino(nupdg)) prob[0]=0.;
  if (pdg::IsAntiNeutrino(nupdg)) prob[2]=0.;

  GHepParticle * target_nucleus = event->TargetNucleus();
  assert(target_nucleus);
  int Z = target_nucleus->Z();
  int A = target_nucleus->A();
  int np = Z;
  int nn = A-Z;

  // select di-nucleon
  bool selected    = false;
  int  cluster_pdg = 0;
  RandomGen * rnd = RandomGen::Instance();
  while(!selected) {
     double prob_gen = rnd->RndGen().Rndm();
     double prob_sum = 0;
     for(int ic = 0; ic < nc; ic++) {
        prob_sum += prob[ic];
        if(prob_gen < prob_sum) {
           cluster_pdg = clusters[ic];
           break;
        }
     }
     // allow selected?
     if ( cluster_pdg==kPdgClusterNN && nn>=2          ) selected=true;
     if ( cluster_pdg==kPdgClusterNP && nn>=1 && np>=1 ) selected=true;
     if ( cluster_pdg==kPdgClusterPP && np>=2          ) selected=true;
  }

  // vtx position in the nucleus
  GHepParticle * neutrino = event->Probe();
  assert(neutrino);
  const TLorentzVector v4(*neutrino->X4());

  // di-nucleon 4-momentum
  double mc = PDGLibrary::Instance()->Find(cluster_pdg)->Mass();
  const TLorentzVector p4(0.,0.,0., mc); 

  // mother particle
  int momidx = event->TargetNucleusPosition();

  // add to the event record
  LOG("MEC", pINFO) 
    << "Adding nucleon cluster [pdgc = " << cluster_pdg << "]";

  event->AddParticle(
     cluster_pdg, kIStNucleonTarget, momidx,-1,-1,-1, p4, v4);
}
//___________________________________________________________________________
void MECGenerator::AddTargetRemnant(GHepRecord * event) const
{
// Add the remnant nucleus (= initial nucleus - nucleon cluster) in the
// event record.

  GHepParticle * target  = event->TargetNucleus();
  GHepParticle * cluster = event->Particle(2);

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

  int momidx = event->TargetNucleusPosition();
  event->AddParticle(ipdgc,kIStStableFinalState, momidx,-1,-1,-1, p4,v4);  
}
//___________________________________________________________________________
void MECGenerator::GenerateFermiMomentum(GHepRecord * event) const
{
// Generate the initial state di-nucleon cluster 4-momentum.
// Draw Fermi momenta for each of the two nucleons.
// Sum the two Fermi momenta to calculate the di-nucleon momentum.
// For simplicity, keep the di-nucleon cluster on the mass shell.
//
  GHepParticle * target_nucleus = event->TargetNucleus();
  assert(target_nucleus);
  GHepParticle * nucleon_cluster = event->Particle(2);
  assert(nucleon_cluster);
  GHepParticle * remnant_nucleus = event->RemnantNucleus();
  assert(remnant_nucleus);
       
  // generate a Fermi momentum for each nucleon

  Target tgt(target_nucleus->Pdg());
  PDGCodeList pdgv = this->NucleonClusterConstituents(nucleon_cluster->Pdg());
  assert(pdgv.size()==2);
  tgt.SetHitNucPdg(pdgv[0]);
  fNuclModel->GenerateNucleon(tgt);
  TVector3 p3a = fNuclModel->Momentum3();
  tgt.SetHitNucPdg(pdgv[1]);
  fNuclModel->GenerateNucleon(tgt);
  TVector3 p3b = fNuclModel->Momentum3();
    
  LOG("FermiMover", pINFO)
     << "1st nucleon (code = " << pdgv[0] << ") generated momentum: ("
     << p3a.Px() << ", " << p3a.Py() << ", " << p3a.Pz() << "), "
     << "|p| = " << p3a.Mag();
  LOG("FermiMover", pINFO)
     << "2nd nucleon (code = " << pdgv[1] << ") generated momentum: ("
     << p3b.Px() << ", " << p3b.Py() << ", " << p3b.Pz() << "), "
     << "|p| = " << p3b.Mag();

  // calcute nucleon cluster momentum

  TVector3 p3 = p3a + p3b;
    
  LOG("FermiMover", pINFO)
     << "di-nucleon cluster momentum: ("
     << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << "), "
     << "|p| = " << p3.Mag();

  // target (initial) nucleus and nucleon cluster mass

  double Mi  = PDGLibrary::Instance()->Find(target_nucleus->Pdg() )-> Mass(); // initial nucleus mass
  double M2n = PDGLibrary::Instance()->Find(nucleon_cluster->Pdg())-> Mass(); // nucleon cluster mass

  // nucleon cluster energy

  double EN = TMath::Sqrt(p3.Mag2() + M2n*M2n);

  // calculate & set the remnant nucleus and nucleon cluster 4-momenta

  TLorentzVector p4nclust   (   p3.Px(),    p3.Py(),    p3.Pz(),  EN   );
  TLorentzVector p4remnant   (-1*p3.Px(), -1*p3.Py(), -1*p3.Pz(), Mi-EN);
       
  nucleon_cluster->SetMomentum(p4nclust);
  remnant_nucleus->SetMomentum(p4remnant);
} 
//___________________________________________________________________________ 
void MECGenerator::SelectKinematics(GHepRecord * event) const
{
// Select interaction kinematics using the rejection method.
//

  // Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  Interaction * interaction = event->Summary();

  // **** NOTE / TODO:
  // **** Hardcode bogus limits for the time-being
  // **** Should be able to get limits via Interaction::KPhaseSpace
  double Q2min =  0.01;
  double Q2max =  8.00;
  double Wmin  =  1.88;
  double Wmax  =  4.50;

  // Scan phase-space for the maximum differential cross section 
  // at the current neutrino energy
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
  RandomGen * rnd = RandomGen::Instance();
  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("MEC", pWARN)
           << "Couldn't select a valid W, Q^2 pair after " 
           << iter << " iterations";
        event->EventFlags()->SetBitNumber(kKineGenErr, true);
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

     // **** NOTE / TODO:
     // **** Forcing a characteristic set of kinematical variables for test purposes.
     // **** Remove once the differential cross section model is implemented
     //     gQ2 = 1.2;
     //     gW  = 2.5;
     //     accept = true;

     // If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("MEC", pINFO) << "Selected: Q^2 = " << gQ2 << ", W = " << gW;

        double Ev = this->EnuAtNucleonClusterRestFrame(event);

        double gx = 0;
        double gy = 0;
        kinematics::WQ2toXY(Ev,2*kNucleonMass,gW,gQ2,gx,gy);

        LOG("MEC", pINFO) << "x = " << gx << ", y = " << gy;

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
void MECGenerator::AddFinalStateLepton(GHepRecord * event) const
{
// **** NOTE / TODO:
// **** Something is broken here --> Odd lepton kinematics.

// Add the final-state primary lepton in the event record.
// Compute its 4-momentum based on the selected interaction kinematics.
//
  Interaction * interaction = event->Summary();
  const InitialState & init_state = interaction->InitState();

  // Look-up selected kinematics
  double Q2 = interaction->Kine().Q2(true);
  double y  = interaction->Kine().y(true);

  // Auxiliary params
  double Ev = this->EnuAtNucleonClusterRestFrame(event);

  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  double El  = (1-y)*Ev;
  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)

  LOG("MEC", pNOTICE)
          << "fsl: E = " << El << ", |p//| = " << plp << ", |pT| = " << plt;

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);

  // Take a unit vector along the neutrino direction
  TVector3 unit_nudir = event->Probe()->P4()->Vect().Unit();

  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the LAB
  TVector3 p3l(pltx,plty,plp);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in LAB
  TLorentzVector p4l(p3l,El);

  // Figure out the final-state primary lepton PDG code
  int pdgc = interaction->FSPrimLepton()->PdgCode();

  // Lepton 4-position (= interacton vtx)
  TLorentzVector v4(*event->Probe()->X4());

  int momidx = event->ProbePosition();
  event->AddParticle(
    pdgc, kIStStableFinalState, momidx, -1, -1, -1, p4l, v4);
}
//___________________________________________________________________________
void MECGenerator::RecoilNucleonCluster(GHepRecord * event) const 
{
  // get di-nucleon cluster 4-momentum
  GHepParticle * nucleon_cluster = event->Particle(2);
  assert(nucleon_cluster);
  TLorentzVector p4cluster(*nucleon_cluster->GetP4());

  // get neutrino 4-momentum
  GHepParticle * neutrino = event->Probe();
  assert(neutrino);
  TLorentzVector p4v(*neutrino->P4());

  // get final state primary lepton 4-momentum
  GHepParticle * fsl = event->FinalStatePrimaryLepton();
  assert(fsl);
  TLorentzVector p4l(*fsl->P4());

  // calculate 4-momentum transfer
  TLorentzVector q = p4v - p4l;

  // calculate recoil nucleon cluster 4-momentum
  TLorentzVector p4cluster_recoil = p4cluster + q;

  // vtx
  TLorentzVector v4(*neutrino->X4());

  // add to the event record
  event->AddParticle(
    nucleon_cluster->Pdg(), kIStDecayedState, 
    2, -1, -1, -1, p4cluster_recoil, v4);
}
//___________________________________________________________________________ 
void MECGenerator::DecayNucleonCluster(GHepRecord * event) const 
{
// Perform a phase-space decay of the nucleon cluster and add its decay
// products in the event record
//
  LOG("MEC", pINFO) << "Decaying nucleon cluster...";

  // get di-nucleon cluster
  int nucleon_cluster_id = 5;
  GHepParticle * nucleon_cluster = event->Particle(nucleon_cluster_id);
  assert(nucleon_cluster);

  // get decay products
  PDGCodeList pdgv = this->NucleonClusterConstituents(nucleon_cluster->Pdg());
  LOG("MEC", pINFO) << "Decay product IDs: " << pdgv;

  // Get the decay product masses
  vector<int>::const_iterator pdg_iter;
  int i = 0;
  double * mass = new double[pdgv.size()];
  double   sum  = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
    int pdgc = *pdg_iter;
    double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
    mass[i++] = m;
    sum += m;
  }

  LOG("MEC", pINFO) 
    << "Performing a phase space decay to "
    << pdgv.size() << " particles / total mass = " << sum;

  TLorentzVector * p4d = nucleon_cluster->GetP4();
  TLorentzVector * v4d = nucleon_cluster->GetX4();

  LOG("MEC", pINFO) 
    << "Decaying system p4 = " << utils::print::P4AsString(p4d);

  // Set the decay
  bool permitted = fPhaseSpaceGenerator.SetDecay(*p4d, pdgv.size(), mass);
  if(!permitted) {
     LOG("MEC", pERROR) 
       << " *** Phase space decay is not permitted \n"
       << " Total particle mass = " << sum << "\n"
       << " Decaying system p4 = " << utils::print::P4AsString(p4d);
     // clean-up 
     delete [] mass;
     delete p4d;
     delete v4d; 
     // throw exception
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Decay not permitted kinematically");
     exception.SwitchOnFastForward();
     throw exception;
  }

  // Get the maximum weight
  double wmax = -1;
  for(int i=0; i<200; i++) {
     double w = fPhaseSpaceGenerator.Generate();   
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  wmax *= 2;

  LOG("MEC", pNOTICE) 
     << "Max phase space gen. weight = " << wmax;

  RandomGen * rnd = RandomGen::Instance();
  bool accept_decay=false;
  unsigned int itry=0;
  while(!accept_decay) 
  {
     itry++;

     if(itry > controls::kMaxUnweightDecayIterations) {
       // report, clean-up and return
       LOG("MEC", pWARN) 
           << "Couldn't generate an unweighted phase space decay after " 
           << itry << " attempts";
       // clean up
       delete [] mass;
       delete p4d;
       delete v4d;
       // throw exception
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select decay after N attempts");
       exception.SwitchOnFastForward();
       throw exception;
     }
     double w  = fPhaseSpaceGenerator.Generate();   
     if(w > wmax) {
        LOG("MEC", pWARN) 
           << "Decay weight = " << w << " > max decay weight = " << wmax;
     }
     double gw = wmax * rnd->RndDec().Rndm();
     accept_decay = (gw<=w);

     LOG("MEC", pINFO) 
        << "Decay weight = " << w << " / R = " << gw 
        << " - accepted: " << accept_decay;

  } //!accept_decay

  // Insert the decay products in the event record
  TLorentzVector v4(*v4d); 
  GHepStatus_t ist = kIStHadronInTheNucleus;
  int idp = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
     int pdgc = *pdg_iter;
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(idp);
     event->AddParticle(pdgc, ist, nucleon_cluster_id,-1,-1,-1, *p4fin, v4);
     idp++;
  }

  // Clean-up
  delete [] mass;
  delete p4d;
  delete v4d;
}
//___________________________________________________________________________
PDGCodeList MECGenerator::NucleonClusterConstituents(int pdgc) const
{
  bool allowdup = true;
  PDGCodeList pdgv(allowdup);

  if(pdgc == kPdgClusterNN) { 
     pdgv.push_back(kPdgNeutron);
     pdgv.push_back(kPdgNeutron);
  }
  else
  if(pdgc == kPdgClusterNP) { 
     pdgv.push_back(kPdgNeutron);
     pdgv.push_back(kPdgProton);
  }
  else
  if(pdgc == kPdgClusterPP) { 
     pdgv.push_back(kPdgProton);
     pdgv.push_back(kPdgProton);
  }
  else 
  {
     LOG("MEC", pERROR) 
        << "Unknown di-nucleon cluster PDG code (" << pdgc << ")";
  }
 
  return pdgv;
}
//___________________________________________________________________________
double MECGenerator::EnuAtNucleonClusterRestFrame(GHepRecord * event) const
{
  GHepParticle * neutrino = event->Particle(0);
  assert(neutrino);
  TLorentzVector p4v(*neutrino->P4());

  GHepParticle * nucleon_cluster = event->Particle(2);
  assert(nucleon_cluster);
  double bx = nucleon_cluster->Px() / nucleon_cluster->Energy();
  double by = nucleon_cluster->Py() / nucleon_cluster->Energy();
  double bz = nucleon_cluster->Pz() / nucleon_cluster->Energy();

  p4v.Boost(-bx,-by,-bz);
  double Ev = p4v.Energy();

  LOG("MEC", pINFO) << "Ev at di-nucleon rest frame = " << Ev;

  return Ev;
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
  fNuclModel = 0;
      
  RgKey nuclkey = "NuclearModel";
  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);
}
//___________________________________________________________________________

