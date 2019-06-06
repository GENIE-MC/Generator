//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab               

         Steve Dytman <dytman+ \at pitt.edu>
         Pittsburgh University

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 22, 2008 - CA
   Skeleton was first added in version 2.5.1
 @ Nov 24-30, 2010 - CA
   Major development leading to the first complete version of the generator.
 @ Nov 20, 2015 - CA, SD  
   Add proper exception handling for failure of phase space decay.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Multinucleon/EventGen/MECGenerator.h"

#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/Multinucleon/XSection/MECHadronTensor.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/PrintUtils.h"

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
  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();
  if (fXSecModel->Id().Name() == "genie::EmpiricalMECPXSec2015") {
      this -> AddTargetRemnant      (event); /// shortly, this will be handled by the InitialStateAppender module
      this -> GenerateFermiMomentum(event);
      this -> SelectEmpiricalKinematics(event);
      // TODO: test removing `AddFinalStateLepton` and replacing it with
      //    PrimaryLeptonGenerator::ProcessEventRecord(evrec);
      this -> AddFinalStateLepton(event);
      // TODO: maybe put `RecoilNucleonCluster` in a `MECHadronicSystemGenerator` class,
      // name it `GenerateEmpiricalDiNucleonCluster` or something...
      this -> RecoilNucleonCluster(event);
      // TODO: `DecayNucleonCluster` should probably be in `MECHadronicSystemGenerator`,
      // if we make that...
      this -> DecayNucleonCluster(event);
  } else if (fXSecModel->Id().Name() == "genie::NievesSimoVacasMECPXSec2016") {
      this -> SelectNSVLeptonKinematics(event);
      this -> AddTargetRemnant(event);
      this -> GenerateNSVInitialHadrons(event);
      // Note: this method in `MECTensor/MECTensorGenerator.cxx` appeared to be a straight
      // copy of an earlier version of the `DecayNucleonCluster` method here - but, watch
      // for this...
      this -> DecayNucleonCluster(event);
  }
  else {
      LOG("MECGenerator",pFATAL) <<
          "ProcessEventRecord >> Cannot calculate kinematics for " <<
          fXSecModel->Id().Name();
  }


}
//___________________________________________________________________________
void MECGenerator::AddTargetRemnant(GHepRecord * event) const
{
// Add the remnant nucleus (= initial nucleus - nucleon cluster) in the
// event record.

  GHepParticle * target  = event->TargetNucleus();
  GHepParticle * cluster = event->HitNucleon();

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

  /*
  if( A == 2 && Z == 2){
    // residual nucleus was just two protons, not a nucleus we know.
    // this might not conserve energy...
    event->AddParticle(kPdgProton,kIStStableFinalState, momidx,-1,-1,-1, p4,v4);  
    // v4,v4 because I'm lazy, give the four momentum to one of the protons, not the other
    event->AddParticle(kPdgProton,kIStStableFinalState, momidx,-1,-1,-1, v4,v4);  
  } else if ( A == 2 && Z == 0){
    // residual nucleus was just two neutrons, not a nucleus we know.
    // this might not conserve energy...
    event->AddParticle(kPdgNeutron,kIStStableFinalState, momidx,-1,-1,-1, p4,v4);  
    // v4,v4 because I'm lazy, give the four momentum to one of the protons, not the other
    event->AddParticle(kPdgNeutron,kIStStableFinalState, momidx,-1,-1,-1, v4,v4);  
  } else {
    // regular nucleus, including deuterium
    event->AddParticle(ipdgc,kIStStableFinalState, momidx,-1,-1,-1, p4,v4);  
  }
  */

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
  GHepParticle * nucleon_cluster = event->HitNucleon();
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

  // set the remnant nucleus and nucleon cluster 4-momenta at GHEP record

  TLorentzVector p4nclust   (   p3.Px(),    p3.Py(),    p3.Pz(),  EN   );
  TLorentzVector p4remnant   (-1*p3.Px(), -1*p3.Py(), -1*p3.Pz(), Mi-EN);
       
  nucleon_cluster->SetMomentum(p4nclust);
  remnant_nucleus->SetMomentum(p4remnant);

  // set the nucleon cluster 4-momentum at the interaction summary 

  event->Summary()->InitStatePtr()->TgtPtr()->SetHitNucP4(p4nclust);
} 
//___________________________________________________________________________ 
void MECGenerator::SelectEmpiricalKinematics(GHepRecord * event) const
{
// Select interaction kinematics using the rejection method.
//

  // Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  Interaction * interaction = event->Summary();
  double Ev = interaction->InitState().ProbeE(kRfHitNucRest);

  // **** NOTE / TODO:
  // **** Hardcode bogus limits for the time-being
  // **** Should be able to get limits via Interaction::KPhaseSpace
  double Q2min =  0.01;
  double Q2max =  8.00;
  double Wmin  =  1.88;
  double Wmax  =  3.00;

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
  LOG("MEC", pNOTICE) << "xsec_max (E = " << Ev << " GeV) = " << xsec_max;

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

     // If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("MEC", pINFO) << "Selected: Q^2 = " << gQ2 << ", W = " << gW;
        double gx = 0;
        double gy = 0;
 	//  More accurate calculation of the mass of the cluster than 2*Mnucl
 	int nucleon_cluster_pdg = interaction->InitState().Tgt().HitNucPdg();
 	double M2n = PDGLibrary::Instance()->Find(nucleon_cluster_pdg)->Mass(); 
 	kinematics::WQ2toXY(Ev,M2n,gW,gQ2,gx,gy);

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
// Add the final-state primary lepton in the event record.
// Compute its 4-momentum based on the selected interaction kinematics.
//
  Interaction * interaction = event->Summary();

  // apapadop: The boost back to the lab frame was missing, that is included now with the introduction of the beta factor
  const InitialState & init_state = interaction->InitState();
  const TLorentzVector & pnuc4 = init_state.Tgt().HitNucP4(); //[@LAB]
  TVector3 beta = pnuc4.BoostVector();

  // apapadop: Boosting the incoming neutrino to the NN-cluster rest frame
  // Neutrino 4p
  TLorentzVector * p4v = event->Probe()->GetP4(); // v 4p @ LAB
  p4v->Boost(-1.*beta);                           // v 4p @ NN-cluster rest frame

  // Look-up selected kinematics
  double Q2 = interaction->Kine().Q2(true);
  double y  = interaction->Kine().y(true);

  // Auxiliary params
  double Ev  = interaction->InitState().ProbeE(kRfHitNucRest);
  LOG("MEC", pNOTICE) << "neutrino energy = " << Ev;
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
  // apapadop: WE NEED THE UNIT VECTOR ALONG THE NEUTRINO DIRECTION IN THE NN-CLUSTER REST FRAME, NOT IN THE LAB FRAME
 TVector3 unit_nudir = p4v->Vect().Unit();      //We use this one, which is in the NN-cluster rest frame
  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the LAB
  TVector3 p3l(pltx,plty,plp);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in LAB
  TLorentzVector p4l(p3l,El);

  // apapadop: Boost final state primary lepton to the lab frame
  p4l.Boost(beta); // active Lorentz transform

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
  // get di-nucleon cluster & its 4-momentum
  GHepParticle * nucleon_cluster = event->HitNucleon();
  assert(nucleon_cluster);
  TLorentzVector * tmp=nucleon_cluster->GetP4();
  TLorentzVector p4cluster(*tmp);
  delete tmp;

  // get neutrino & its 4-momentum
  GHepParticle * neutrino = event->Probe();
  assert(neutrino);
  TLorentzVector p4v(*neutrino->P4());

  // get final state primary lepton & its 4-momentum
  GHepParticle * fsl = event->FinalStatePrimaryLepton();
  assert(fsl);
  TLorentzVector p4l(*fsl->P4());

  // calculate 4-momentum transfer
  TLorentzVector q = p4v - p4l;

  // calculate recoil nucleon cluster 4-momentum
  TLorentzVector p4cluster_recoil = p4cluster + q;

  // work-out recoil nucleon cluster code
  LOG("MEC", pINFO) << "Interaction summary";
  LOG("MEC", pINFO) << *event->Summary();
  int recoil_nucleon_cluster_pdg = event->Summary()->RecoilNucleonPdg();

  // vtx position in nucleus coord system
  TLorentzVector v4(*neutrino->X4());

  // add to the event record
  event->AddParticle(
    recoil_nucleon_cluster_pdg, kIStDecayedState, 
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
     event->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Decay not permitted kinematically");
     exception.SwitchOnStepBack();
     exception.SetReturnStep(0);
     throw exception;
  }

  // Get the maximum weight
  double wmax = -1;
  for(int idec=0; idec<200; idec++) {
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
       event->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select decay after N attempts");
       exception.SwitchOnStepBack();
       exception.SetReturnStep(0);
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
void MECGenerator::SelectNSVLeptonKinematics (GHepRecord * event) const
{
  // -- implementation -- //
  // The IFIC Valencia model can provide three different hadron tensors.
  // The user probably wants all CC QE-like 2p2h events
  // But could in principle get the no-delta component if they want (deactivated incode)
  int FullDeltaNodelta = 1;  // 1:  full, 2:  only delta, 3:  zero delta

  // -- limit the maximum XS for the accept/reject loop -- //
  // 
  // MaxXSec parameters.  This whole calculation could be in it's own function?
  // these need to lead to a number that is safely large enough, or crash the run.
  double XSecMaxPar1 = 2.2504;
  double XSecMaxPar2 = 9.41158;

  // -- Event Properties -----------------------------//
  Interaction * interaction = event->Summary();
  Kinematics * kinematics = interaction->KinePtr();

  double Enu = interaction->InitState().ProbeE(kRfHitNucRest);

  int NuPDG = interaction->InitState().ProbePdg();
  int TgtPDG = interaction->InitState().TgtPdg();
  // interacton vtx
  TLorentzVector v4(*event->Probe()->X4());
  TLorentzVector tempp4(0.,0.,0.,0.);

  // -- Lepton Kinematic Limits ----------------------------------------- //
  double Costh = 0.0; // lepton angle
  double CosthMax = 1.0;
  double CosthMin = -1.0;

  double T = 0.0;  // lepton kinetic energy
  double TMax = std::numeric_limits<double>::max();
  double TMin = 0.0;

  double Plep = 0.0; // lepton 3 momentum
  double Elep = 0.0; // lepton energy
  double LepMass = interaction->FSPrimLepton()->Mass();
  
  double Q0 = 0.0; // energy component of q four vector
  double Q3 = 0.0; // magnitude of transfered 3 momentum
  double Q2 = 0.0; // properly Q^2 (Q squared) - transfered 4 momentum.

  // Set lepton KE TMax for for throwing rndm in the accept/reject loop.
  // We can accidentally set it too high, because the xsec will return zero.
  // This way if someone reuses this code, they are not tripped up by it.
  TMax = Enu - LepMass;

  // Set Tmin for throwing rndm in the accept/reject loop
  // the hadron tensors we expect will be limited in q3
  // therefore also the outgoing lepton KE can't be too low or costheta too backward
  // make the accept/reject loop more efficient by using Min values.
  if(Enu < fQ3Max){
    TMin = 0 ;
    CosthMin = -1 ; 
  } else {
    TMin = TMath::Sqrt(TMath::Power(LepMass, 2) + TMath::Power((Enu - fQ3Max), 2)) - LepMass;
    CosthMin = TMath::Sqrt(1 - TMath::Power((fQ3Max / Enu ), 2));
  }

  // The accept/reject loop tests a rand against a maxxsec - must scale with A.
  int NuclearA = 12;
  double NuclearAfactorXSecMax = 1.0;
  if (TgtPDG != kPdgTgtC12) {
    if (TgtPDG > kPdgTgtFreeN && TgtPDG) {
      NuclearA = pdg::IonPdgCodeToA(TgtPDG);
      // The QE-like portion scales as A, but the Delta portion increases faster, not simple.
      // so this gives additional safety factor.  Remember, we need a safe max, not precise max.
      if (NuclearA < 12) NuclearAfactorXSecMax *= NuclearA / 12.0;
      else NuclearAfactorXSecMax *= TMath::Power(NuclearA/12.0, 1.4);
    } 
    else {
      LOG("MEC", pERROR) << "Trying to scale XSecMax for larger nuclei, but "
          << TgtPDG << " isn't a nucleus?";
      assert(false);
    }
  }
  
  // -- Generate and Test the Kinematics----------------------------------//

  RandomGen * rnd = RandomGen::Instance();
  bool accept = false;
  unsigned int iter = 0;

  // loop over different (randomly) selected T and Costh
  while (!accept) {
      iter++;
      if(iter > kRjMaxIterations) {
          // error if try too many times
          LOG("MEC", pWARN)
              << "Couldn't select a valid Tmu, CosTheta pair after " 
              << iter << " iterations";
          event->EventFlags()->SetBitNumber(kKineGenErr, true);
          genie::exceptions::EVGThreadException exception;
          exception.SetReason("Couldn't select lepton kinematics");
          exception.SwitchOnFastForward();
          throw exception;
      }

      // generate random kinetic energy T and Costh
      T = TMin + (TMax-TMin)*rnd->RndKine().Rndm();
      Costh = CosthMin + (CosthMax-CosthMin)*rnd->RndKine().Rndm();

      // Calculate useful values for judging this choice
      Plep = TMath::Sqrt( T * (T + (2.0 * LepMass)));  // ok is sqrt(E2 - m2)
      Q3 = TMath::Sqrt(Plep*Plep + Enu*Enu - 2.0 * Plep * Enu * Costh);

      // Don't bother doing hard work if the selected Q3 is greater than Q3Max
      if (Q3 < fQ3Max){

          kinematics->SetKV(kKVTl, T);
          kinematics->SetKV(kKVctl, Costh);

          // decide whether to accept or reject these kinematics
          // AND set the chosen two-nucleon system

          // to save time, use a pre-calculated max cross-section XSecMax
          // it doesn't matter what it is, as long as it is big enough.
          // RIK asks can XSecMax can be pushed to the q0q3 part of the calculation
          // where the XS doesn't depend much on Enu.
          // instead, this implementation uses a rough dependence on log10(Enu).
          // starting around 0.5 GeV, the log10(max) is linear vs. log10(Enu)
          // 1.10*TMath::Power(10.0, 2.2504 * TMath::Log10(Enu) - 9.41158)

          if (FullDeltaNodelta == 1){ 
              // this block for the user who wants all CC QE-like 2p2h events

              // extract xsecmax from the spline making process for C12 and other nuclei.
              //  plot Log10(E) on horizontal and Log10(xsecmax) vertical
              //  and fit a line.  Use that plus 1.35 safety factors to limit the accept/reject loop.
              double XSecMax = 1.35 * TMath::Power(10.0, XSecMaxPar1 * TMath::Log10(Enu) - XSecMaxPar2);
              if (NuclearA > 12) XSecMax *=  NuclearAfactorXSecMax;  // Scale it by A, precomputed above.

              LOG("MEC", pDEBUG) << " T, Costh: " << T << ", " << Costh ;


              // We need four different cross sections. Right now, pursue the
              // inelegant method of calling XSec four times - there is
              // definitely some runtime inefficiency here, but it is not awful

              // first, get delta-less all
              if (NuPDG > 0) {
                  interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgClusterNN);
              }
              else {
                  interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgClusterPP);
              }
              double XSec = fXSecModel->XSec(interaction, kPSTlctl);
              // now get all with delta
              interaction->ExclTagPtr()->SetResonance(genie::kP33_1232);
              double XSecDelta = fXSecModel->XSec(interaction, kPSTlctl);
              // get PN with delta
              interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgClusterNP);
              double XSecDeltaPN = fXSecModel->XSec(interaction, kPSTlctl);
              // now get delta-less PN
              interaction->ExclTagPtr()->SetResonance(genie::kNoResonance);
              double XSecPN = fXSecModel->XSec(interaction, kPSTlctl);

              if (XSec > XSecMax) {
                  LOG("MEC", pERROR) << "XSec is > XSecMax for nucleus " << TgtPDG << " " 
				   << XSec << " > " << XSecMax 
				   << " don't let this happen.";
              }
              assert(XSec <= XSecMax);
              accept = XSec > XSecMax*rnd->RndKine().Rndm();
              LOG("MEC", pINFO) << "Xsec, Max, Accept: " << XSec << ", " 
                  << XSecMax << ", " << accept; 

              if(accept){
                  // If it passes the All cross section we still need to do two things:
                  // * Was the initial state pn or not?
                  // * Do we assign the reaction to have had a Delta on the inside?

                  // PDD means from the part of the XSec with an internal Delta line
                  // that (at the diagram level) did not produce a pion in the final state.

                  bool isPDD = false;

                  // Find out if we should use a pn initial state
                  double myrand = rnd->RndKine().Rndm();
                  double pnFraction = XSecPN / XSec;
                  LOG("MEC", pDEBUG) << "Test for pn: xsec_pn = " << XSecPN 
                      << "; xsec = " << XSec 
                      << "; pn_fraction = " << pnFraction
                      << "; random number val = " << myrand;

                  if (myrand <= pnFraction) {
                      // yes it is, add a PN initial state to event record
                      event->AddParticle(kPdgClusterNP, kIStNucleonTarget,
                              1, -1, -1, -1, tempp4, v4);
                      interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgClusterNP);

                      // Its a pn, so test for Delta by comparing DeltaPN/PN
                      if (rnd->RndKine().Rndm() <= XSecDeltaPN / XSecPN) {
                          isPDD = true;
                      }
                  }
                  else {
                      // no it is not a PN, add either NN or PP initial state to event record.
                      if (NuPDG > 0) {
                          event->AddParticle(kPdgClusterNN, kIStNucleonTarget,
                                  1, -1, -1, -1, tempp4, v4);
                          interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgClusterNN);
                      }
                      else {
                          event->AddParticle(kPdgClusterPP, kIStNucleonTarget,
                                  1, -1, -1, -1, tempp4, v4);
                          interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgClusterPP); 
                      }
                      // its not pn, so test for Delta (XSecDelta-XSecDeltaPN)/(XSec-XSecPN)
                      // right, both numerator and denominator are total not pn.
                      if (rnd->RndKine().Rndm() <=
                              (XSecDelta - XSecDeltaPN) / (XSec - XSecPN)) {
                          isPDD = true;
                      }
                  }

                  // now test whether we tagged this as a pion event
                  // and assign that fact to the Exclusive State tag
                  // later, we can query const XclsTag & xcls = interaction->ExclTag() 
                  if (isPDD){
                      interaction->ExclTagPtr()->SetResonance(genie::kP33_1232);
                  }


              } // end if accept
          } // end if delta == 1

          /* One can make simpler versions of the above for the
             FullDeltaNodelta == 2 (only delta)
             or
             FullDeltaNodelta == 3 (set Delta FF = 1, lose interference effect).
             but I (Rik) don't see what the use-case is for these, genratorly speaking.
             */

      }// end if passes q3 test
  } // end while

  // -- finish lepton kinematics
  // If the code got here, then we accepted some kinematics
  // and we can proceed to generate the final state.

  // define coordinate system wrt neutrino: z along neutrino, xy perp

  // Cos theta gives us z, the rest in xy:
  double PlepZ = Plep * Costh;
  double PlepXY = Plep * TMath::Sqrt(1. - TMath::Power(Costh,2));

  // random rotation about unit vector for phi direction
  double phi= 2 * kPi * rnd->RndLep().Rndm();
  // now fill x and y from PlepXY
  double PlepX = PlepXY * TMath::Cos(phi);
  double PlepY = PlepXY * TMath::Sin(phi);

  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the LAB
  TVector3 unit_nudir = event->Probe()->P4()->Vect().Unit();
  TVector3 p3l(PlepX, PlepY, PlepZ);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in LAB
  Elep = TMath::Sqrt(LepMass*LepMass + PlepX*PlepX + PlepY*PlepY + PlepZ*PlepZ); 
  TLorentzVector p4l(p3l,Elep);

  // Figure out the final-state primary lepton PDG code
  int pdgc = interaction->FSPrimLepton()->PdgCode();
  int momidx = event->ProbePosition();

  // -- Store Values ------------------------------------------//
  // -- Interaction: Q2
  Q0 = Enu - Elep;
  Q2 = Q3*Q3 - Q0*Q0;
  double gy = Q0 / Enu;
  double gx = kinematics::Q2YtoX(Enu, 2 * kNucleonMass, Q2, gy);
  double gW = kinematics::XYtoW(Enu, 2 * kNucleonMass, gx, gy);

  interaction->KinePtr()->SetQ2(Q2, true);
  interaction->KinePtr()->Sety(gy, true);
  interaction->KinePtr()->Setx(gx, true);
  interaction->KinePtr()->SetW(gW, true);
  interaction->KinePtr()->SetFSLeptonP4(p4l);
  // in later methods
  // will also set the four-momentum and W^2 of the hadron system.

  // -- Lepton
  event->AddParticle( pdgc, kIStStableFinalState, momidx, -1, -1, -1, p4l, v4);

  LOG("MEC",pDEBUG) << "~~~ LEPTON DONE ~~~";
}
//___________________________________________________________________________
void MECGenerator::GenerateNSVInitialHadrons(GHepRecord * event) const
{
    // We need a kinematic limits accept/reject loop here, so generating the
    // initial hadrons is combined with generating the recoil hadrons...

    LOG("MEC",pDEBUG) << "Generate Initial Hadrons - Start";

    Interaction * interaction = event->Summary();
    GHepParticle * neutrino = event->Probe();
    assert(neutrino);
    TLorentzVector p4nu(*neutrino->P4());

    // get final state primary lepton & its 4-momentum
    GHepParticle * fsl = event->FinalStatePrimaryLepton();
    assert(fsl);
    TLorentzVector p4l(*fsl->P4());

    // calculate 4-momentum transfer from these
    TLorentzVector Q4 = p4nu - p4l;

    // get the target two-nucleon cluster and nucleus.
    // the remnant nucleus is apparently set, except for its momentum.
    GHepParticle * target_nucleus = event->TargetNucleus();
    assert(target_nucleus);
    GHepParticle * initial_nucleon_cluster = event->HitNucleon();
    assert(initial_nucleon_cluster);
    GHepParticle * remnant_nucleus = event->RemnantNucleus();
    assert(remnant_nucleus);

    // -- make a two-nucleon system, then give them some momenta.

    // instantiate an empty local target nucleus, so I can use existing methods
    // to get a momentum from the prevailing Fermi-motion distribution 
    Target tgt(target_nucleus->Pdg());

    // NucleonClusterConstituents is an implementation within this class, called with this
    // It using the nucleon cluster from the earlier tests for a pn state,
    // the method returns a vector of pdgs, which hopefully will be of size two.

    PDGCodeList pdgv = this->NucleonClusterConstituents(initial_nucleon_cluster->Pdg());
    assert(pdgv.size()==2);

    // These things need to be saved through to the end of the accept loop.
    bool accept = false;
    TVector3 p31i;
    TVector3 p32i;
    unsigned int iter = 0;

    int initial_nucleon_cluster_pdg = initial_nucleon_cluster->Pdg();
    int final_nucleon_cluster_pdg = 0;
    if (neutrino->Pdg() > 0) {
        if (initial_nucleon_cluster->Pdg() == kPdgClusterNP) {
            final_nucleon_cluster_pdg = kPdgClusterPP;
        }
        else if (initial_nucleon_cluster->Pdg() == kPdgClusterNN) {
            final_nucleon_cluster_pdg = kPdgClusterNP;
        }
        else {
            LOG("MEC", pERROR) << "Wrong pdg for a CC neutrino MEC interaction" 
                << initial_nucleon_cluster->Pdg();
        }
    } 
    else if (neutrino->Pdg() < 0) {
        if (initial_nucleon_cluster->Pdg() == kPdgClusterNP) {
            final_nucleon_cluster_pdg = kPdgClusterNN;
        }
        else if (initial_nucleon_cluster->Pdg() == kPdgClusterPP) {
            final_nucleon_cluster_pdg = kPdgClusterNP;
        }
        else {
            LOG("MEC", pERROR) << "Wrong pdg for a CC anti-neutrino MEC interaction" 
                << initial_nucleon_cluster->Pdg();
        }
    }

    TLorentzVector p4initial_cluster;
    TLorentzVector p4final_cluster;
    TLorentzVector p4remnant_nucleus;
    double removalenergy1;
    double removalenergy2;

    //===========================================================================
    // Choose two nucleons from the prevailing fermi-motion distribution.
    // Some produce kinematically unallowed combinations initial cluster and Q2
    // Find out, and if so choose them again with this accept/reject loop.
    // Some kinematics are especially tough 
    while(!accept){
        iter++;
        if(iter > kRjMaxIterations*1000) {
            // error if try too many times
            LOG("MEC", pWARN)
                << "Couldn't select a valid W, Q^2 pair after " 
                << iter << " iterations";
            event->EventFlags()->SetBitNumber(kKineGenErr, true);
            genie::exceptions::EVGThreadException exception;
            exception.SetReason("Couldn't select initial hadron kinematics");
            exception.SwitchOnFastForward();
            throw exception;
        }

        // generate nucleons
        // tgt is a Target object for local use, just waiting to be filled.
        // this sets the pdg of each nucleon and its momentum from a Fermi-gas
        // Nieves et al. would use a local Fermi gas here, not this, but ok.
        // so momentum from global Fermi gas, local Fermi gas, or spectral function
        // and removal energy ~0.025 GeV, correlated with density, or from SF distribution
        tgt.SetHitNucPdg(pdgv[0]);
        fNuclModel->GenerateNucleon(tgt);
        p31i = fNuclModel->Momentum3();
        removalenergy1 = fNuclModel->RemovalEnergy();
        tgt.SetHitNucPdg(pdgv[1]);
        fNuclModel->GenerateNucleon(tgt);
        p32i = fNuclModel->Momentum3();
        removalenergy2 = fNuclModel->RemovalEnergy();

        // not sure -- could give option to use Nieves q-value here.

        // Now write down the initial cluster four-vector for this choice
        TVector3 p3i = p31i + p32i;
        double mass2 = PDGLibrary::Instance()->Find( initial_nucleon_cluster_pdg )->Mass();
        mass2 *= mass2;
        double energy = TMath::Sqrt(p3i.Mag2() + mass2);
        p4initial_cluster.SetPxPyPzE(p3i.Px(),p3i.Py(),p3i.Pz(),energy);

        // cast the removal energy as the energy component of a 4-vector for later.
        TLorentzVector tLVebind(0., 0., 0., -1.0 * (removalenergy1 + removalenergy2));

        // RIK: You might ask why is this the right place to subtract ebind?
        // It is okay. Physically, I'm subtracting it from q. The energy
        // transfer to the nucleon is 50 MeV less. the energy cost to put this
        // cluster on-shell. What Jan says he does in PRC.86.015504 is this:
        //   The nucleons are assumed to be in a potential well
        //     V = Efermi + 8 MeV.
        //   The Fermi energy is subtracted from each initial-state nucleon
        //   (I guess he does it dynamically based on Ef = p2/2M or so which 
        //   is what we are doing above, on average). Then after the reaction,
        // another 8 MeV is subtracted at that point; a small adjustment to the
        // momentum is needed to keep the nucleon on shell.

        // calculate recoil nucleon cluster 4-momentum (tLVebind is negative)
        p4final_cluster = p4initial_cluster + Q4 + tLVebind;

        // Test if the resulting four-vector corresponds to a high-enough invariant mass.
        // Fail the accept if we couldn't put this thing on-shell.
        if (p4final_cluster.M() < 
                PDGLibrary::Instance()->Find(final_nucleon_cluster_pdg )->Mass()) {
            accept = false;
        } else {
            accept = true;
        }

    }  // end accept loop

    // we got here if we accepted the final state two-nucleon system
    // so now we need to write everything to ghep

    // First the initial state nucleons.
    initial_nucleon_cluster->SetMomentum(p4initial_cluster);

    // and the remnant nucleus
    double Mi  = PDGLibrary::Instance()->Find(target_nucleus->Pdg() )-> Mass();
    remnant_nucleus->SetMomentum(-1.0*p4initial_cluster.Px(),
            -1.0*p4initial_cluster.Py(),
            -1.0*p4initial_cluster.Pz(),
            Mi - p4initial_cluster.E() + removalenergy1 + removalenergy2);

    // Now the final nucleon cluster.

    // Getting this v4 is a position, i.e. a position within the nucleus (?)
    // possibly it takes on a non-trivial value only for local Fermi gas
    // or for sophisticated treatments of intranuclear rescattering.
    TLorentzVector v4(*neutrino->X4());

    // Now write the (undecayed) final two-nucleon system
    GHepParticle p1(final_nucleon_cluster_pdg, kIStDecayedState,
            2, -1, -1, -1, p4final_cluster, v4);
    
    //p1.SetRemovalEnergy(removalenergy1 + removalenergy2);
    // The "bound particle" concept applies only to p or n.
    // Instead, add this directly to the remnant nucleon a few lines above.

    // actually, this is not an status1 particle, so it is not picked up
    // by the aggregator. anyway, the aggregator does not run until the very end.

    event->AddParticle(p1);

    interaction->KinePtr()->SetHadSystP4(p4final_cluster);
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

    GetParam( "NSV-Q3Max", fQ3Max ) ;
}
//___________________________________________________________________________

