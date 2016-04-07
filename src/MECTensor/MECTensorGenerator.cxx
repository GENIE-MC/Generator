//________________________________________________________
/*

the accept-reject loop that pulls lepton kinematics and
the pn fraction from the hadron tensor.
Going beyond the model calculation, this code then
assigns the energy and momentum transfer to two nucleons.
the nucleon final state strategy is similar to 
Dytman's MEC code in GENIE and the NuWRO strategy.

2014/09/15 first working version~ J.Schwehr
2015/01/14 all neutrino flavors, oxygen and carbon ~J.Schwehr
2015/02/06 fixed hadron system, added spline generation code ~R. Gran
2016/02/20 general clean-up and fix a problem with PDD ~R. Gran
2016/03/01 modified XSec structure to handle non-isoscalar ~R. Gran

*/
//_________________________________________________________

// -- includes -- //

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
#include "MECTensor/MECTensorGenerator.h"
#include "Numerical/RandomGen.h"
#include "Numerical/BLI2D.h"
#include "Nuclear/NuclearModelI.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"
#include "Utils/PrintUtils.h"
#include "MECTensor/MECLoadHadTensor.h"
#include <iostream>
#include <limits>

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;


//___________________________________________________________________________
MECTensorGenerator::MECTensorGenerator() :
EventRecordVisitorI("genie::MECTensorGenerator")
{
  LOG("MEC", pDEBUG) << "MECTensorGenerator Init";

}
//___________________________________________________________________________
MECTensorGenerator::MECTensorGenerator(string config) :
EventRecordVisitorI("genie::MECTensorGenerator", config)
{
  LOG("MEC", pDEBUG) << "MECTensorGenerator Init";
}
//___________________________________________________________________________
MECTensorGenerator::~MECTensorGenerator()
{

}
//___________________________________________________________________________
void MECTensorGenerator::ProcessEventRecord(GHepRecord * event) const
{

  this -> SelectLeptonKinematics (event);
  this -> AddTargetRemnant       (event);
  this -> GenerateInitialHadrons (event);
  //this -> RecoilNucleonCluster   (event);
  this -> DecayNucleonCluster    (event);

}
//___________________________________________________________________________
void MECTensorGenerator::SelectLeptonKinematics (GHepRecord * event) const
{
  // -- Constants --------------------------------- //
  // these need to be moved or replaced with the "correct" variable

  // -- model -- //
  // these need a configuration file
  // Q3Max needs to be <= the size of the input HadTensor.
  // When we have more than one option available, this should go with the input
  double Q3Max = 1.2;  
  
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

  // -- genie -- //
  // JS: can't find a better way to do this.
  double arbitraryLargePDG= 1002000000; 

  // -- Event Properties -----------------------------//
  Interaction * interaction = event->Summary();
  InitialState * init_state = interaction->InitStatePtr();
  double Enu = interaction->InitState().ProbeE(kRfHitNucRest);
  double LepMass = interaction->FSPrimLepton()->Mass();
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
  
  double Q0 = 0.0; // energy component of q four vector
  double Q3 = 0.0; // magnitude of transfered 3 momentum
  double Q2 = 0.0; // properly Q^2 (Q squared) - transfered 4 momentum.


  // -- load xsec tables --- //
  LOG("MEC", pDEBUG) << "About to run MECLoadHadTensor::Instance from MECTensorGenerator.cxx ";

  // TODO
  // Make a qvalue lookup function to return qvale based on A, nu/nubar
  // And encode it from Juan's table.

  
  MECLoadHadTensor * hadtensor = MECLoadHadTensor::Instance();

  // Set lepton KE TMax for for throwing rndm in the accept/reject loop.
  // We can accidentally set it too high, because the xsec will return zero.
  // This way if someone reuses this code, they are not tripped up by it.
  //double qvalue = hadtensor->Qvalue(TgtPDG,NuPDG); //(0.016827;
  //TMax = Enu - LepMass - qvalue;
  TMax = Enu - LepMass;

  // Set Tmin for throwing rndm in the accept/reject loop
  // the hadron tensors we expect will be limited in q3
  // therefore also the outgoing lepton KE can't be too low or costheta too backward
  // make the accept/reject loop more efficient by using Min values.
  if(Enu < Q3Max){
    TMin = 0 ;
    CosthMin = -1 ; 
  } else {
    TMin = TMath::Sqrt( TMath::Power(LepMass,2) + TMath::Power((Enu-Q3Max),2) ) - LepMass ;
    CosthMin = TMath::Sqrt( 1 - TMath::Power(( Q3Max / Enu ),2) ) ;
  }

  // The accept/reject loop tests a rand against a maxxsec .
  // We should scale it with A.
  // Here I already can know what A is
  // get it and precompute a few things for that later test.

  int NuclearA = 12;
  int NuclearAfactorXSecMax = 1.0;
  if(TgtPDG != kPdgTgtC12){
    if(TgtPDG > kPdgTgtFreeN && TgtPDG < arbitraryLargePDG){
      NuclearA = pdg::IonPdgCodeToA(TgtPDG);
      // The QE-like portion scales as A, but the Delta portion increases faster, not simple.
      // so this gives additional safety factor.  Remember, we need a safe max, not precise max.
      if (NuclearA < 12) NuclearAfactorXSecMax *= NuclearA / 12.0;
      else NuclearAfactorXSecMax *= TMath::Power(NuclearA/12.0, 1.2);
    } else {
      LOG("MEC", pERROR) << "Trying to scale XSecMax for larger nuclei, but " << TgtPDG << " isn't a nucleus?";
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
    if (Q3 < Q3Max){

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
	
	  //double XSecMaxPar1 = 2.2504;  double XSecMaxPar2 = 9.41158;
	  // extract xsecmax from the spline making process for C12 and other nuclei.
	  //  plot Log10(E) on horizontal and Log10(xsecmax) vertical
	  //  and fit a line.  Use that plus 1.25 safety factors to limit the accept/reject loop.
	  double XSecMax = 1.25*TMath::Power(10.0, XSecMaxPar1 * TMath::Log10(Enu) - XSecMaxPar2);
	  XSecMax *=  NuclearAfactorXSecMax;  // Scale it by A, precomputed above.
	
	  LOG("MEC", pDEBUG) << " T, Costh: " << T << ", " << Costh ;

	  double XSecFour[4];

	  hadtensor->XSecFour(TgtPDG, NuPDG, Enu, T, Costh, XSecFour, true);
	  
          //double XSec = hadtensor->XSecFullAll(TgtPDG, NuPDG, Enu, T, Costh);
	  double XSec = XSecFour[0];
          if(XSec > XSecMax) LOG("MEC", pERROR) << "XSec is > XSecMax for nucleus " << TgtPDG << " don't let this happen.";
          assert(XSec <= XSecMax);
          accept = XSec > XSecMax*rnd->RndKine().Rndm();
          LOG("MEC", pINFO) << "Xsec, Max, Accept: " << XSec << ", " << XSecMax << ", " << accept; 


	  if(accept){
	      // If it passes the All cross section we still need to do two things
	      // Was the initial state pn or not?
	      // Do we assign the reaction to have had a Delta on the inside?

	      // This uses the second hadtensor, delta-only to get a second subset Xsec
	      // Set PDD "pionless delta decay" if it passes this also.
	      // PDD means from the part of the XSec with an internal Delta line
	      // that (at the diagram level) did not produce a pion in the final state.
	      // Precompute the two cross sections we need, but throw the random later.
	      //double XSecDelta = hadtensor->XSecDeltaAll(TgtPDG, NuPDG, Enu, T, Costh);
	      //double XSecDeltaPN = hadtensor->XSecDeltapn(TgtPDG, NuPDG, Enu, T, Costh);
	      double XSecDelta = XSecFour[2];
	      double XSecDeltaPN = XSecFour[3];
	      bool isPDD = 0;  // this flag is actually set inside the if blocks.

	      // Find out if we should use a pn initial state
              //double XSecPN= hadtensor->XSecFullpn(TgtPDG, NuPDG, Enu, T, Costh);
	      double XSecPN = XSecFour[1];
	      
              double myrand = rnd->RndKine().Rndm();
              double pnFraction = XSecPN / XSec;
              LOG("MEC", pINFO)<< "RIK test for pn " << XSecPN << " " << XSec << " " << pnFraction << " " << myrand;

              if ( myrand <= pnFraction){
                  // yes it is, add a PN initial state to event record
                  event->AddParticle(kPdgClusterNP, kIStNucleonTarget, 1, -1, -1, -1, tempp4, v4);
                  init_state->TgtPtr()->SetHitNucPdg(kPdgClusterNP);

		  // Its a pn, so test for Delta by comparing DeltaPN/PN
		  if(rnd->RndKine().Rndm() <= XSecDeltaPN/XSecPN)isPDD=1;
		  
              }
              else {
		  // no it is not, add either NN or PP initial state to event record.
                  if (NuPDG > 0) {
                      event->AddParticle(kPdgClusterNN, kIStNucleonTarget, 1, -1, -1, -1, tempp4, v4);
                      init_state->TgtPtr()->SetHitNucPdg(kPdgClusterNN);
                  }
                  else {
                      event->AddParticle(kPdgClusterPP, kIStNucleonTarget, 1, -1, -1, -1, tempp4, v4);
                      init_state->TgtPtr()->SetHitNucPdg(kPdgClusterPP); 
                  }

		  // its not pn, so test for Delta (XSecDelta-XSecDeltaPN)/(XSec-XSecPN)
		  // right, both numerator and denominator are total not pn.
		  if(rnd->RndKine().Rndm() <= (XSecDelta-XSecDeltaPN)/(XSec-XSecPN))isPDD=1;

		  
              }

	      // now test whether we tagged this as a pion event
	      // and assign that fact to the Exclusive State tag
	      // later, we can query const XclsTag & xcls = interaction->ExclTag() 
              if (isPDD){
		XclsTag * xcls = interaction->ExclTagPtr();
		xcls->SetResonance(kP33_1232);		
              }
	      
	  } // end if accept
      }// end if delta ==1
      
      /* One can make simpler versions of the above for the
         FullDeltaNodelta == 2 (only delta)
         or
         FullDeltaNodelta == 3 (set Delta FF = 1, lose interference effect).
         but I don't see what the use-case is for these, genratorly speaking.
      */


    }// end if passes q3 test
  }// end while

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

  // Lepton 4-position (= interacton vtx)
  //TLorentzVector v4(*event->Probe()->X4()); -- moved to start of code

  int momidx = event->ProbePosition();

  // -- SANITY CHECK -- //
  if (Elep > Enu) {
    LOG("MEC",pERROR) << "Energy Problem! Elep " << Elep << "  > Enu " << Enu;
  }
 

  // -- Store Values ------------------------------------------//

  // -- Interaction: Q2
  Q0 = Enu - Elep;
  Q2 = Q3*Q3 - Q0*Q0;

  interaction->KinePtr()->SetQ2(Q2, true);
  interaction->KinePtr()->Sety(Q0/Enu, true);
  interaction->KinePtr()->SetFSLeptonP4(p4l);
  // in later methods
  // will also set the four-momentum and W^2 of the hadron system.

  // -- Lepton
  event->AddParticle( pdgc, kIStStableFinalState, momidx, -1, -1, -1, p4l, v4);
  
  LOG("MEC",pDEBUG) << "~~~ LEPTON DONE ~~~";

}
//____________________________________________________________________

//___________________________________________________________________________
// coppied verbatum //
void MECTensorGenerator::AddTargetRemnant(GHepRecord * event) const
{
// Add the remnant nucleus (= initial nucleus - nucleon cluster) in the
// event record.


  LOG("MEC", pDEBUG) << "Adding Remnant";
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
  event->AddParticle(ipdgc,kIStStableFinalState, momidx,-1,-1,-1, p4,v4);  

  LOG("MEC",pDEBUG) << "Remnant Done";

}
//_________________________________________________________________________
void MECTensorGenerator::GenerateInitialHadrons  (GHepRecord * event) const
{
  // Steve's MEC code separated the GenerateInitialHadrons from the recoil hadrons
  // But we need a kinematic limits accept/reject loop here, so its one method.

  LOG("MEC",pDEBUG) << "Generate Initial Hadrons - Start";

  // -- Inputs: Q4, Fermi Momentum ------------------------N--
  // get neutrino & its 4-momentum

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
  // Get the outgoing cluster mass.
  // the method in the next line does not work for some reason.
  //final_nucleon_cluster_pdg = event->Summary()->RecoilNucleonPdg();
  //std::cout << "RIK final_nucleon " << final_nucleon_cluster_pdg << " and initial " << initial_nucleon_cluster->Pdg() << std::endl;
  //std::cout << "RIK final_nucleon " << final_nucleon_cluster_pdg << " and initial " << initial_nucleon_cluster_pdg << std::endl;
  //heisenbug here.  May be related to shoehorning this into just the 200 or 202 initial states, not 201.

  // get ready to set the pdg code for the final two-nucleon cluster
  // This is lingering here because I once had a problem with RecoilNucleonPdg() method.
  // override that method for now, but really go and fix it.
  if(neutrino->Pdg() > 0){
    if(initial_nucleon_cluster->Pdg() == kPdgClusterNP)final_nucleon_cluster_pdg = kPdgClusterPP;
    else if(initial_nucleon_cluster->Pdg() == kPdgClusterNN)final_nucleon_cluster_pdg = kPdgClusterNP;
    else LOG("MEC", pERROR) << "Wrong pdg for a CC neutrino MEC interaction" << initial_nucleon_cluster->Pdg();
  } else if(neutrino->Pdg() < 0){
    if(initial_nucleon_cluster->Pdg() == kPdgClusterNP)final_nucleon_cluster_pdg = kPdgClusterNN;
    else if(initial_nucleon_cluster->Pdg() == kPdgClusterPP)final_nucleon_cluster_pdg = kPdgClusterNP;
    else LOG("MEC", pERROR) << "Wrong pdg for a CC anti-neutrino MEC interaction" << initial_nucleon_cluster->Pdg();
  }
  
  // sanity check.  might not be needed anymore.
  if(final_nucleon_cluster_pdg != kPdgClusterPP && final_nucleon_cluster_pdg != kPdgClusterNN
     && final_nucleon_cluster_pdg != kPdgClusterNP && final_nucleon_cluster_pdg != initial_nucleon_cluster_pdg){
    LOG("MEC", pERROR) << "Wrong pdg for a CC neutrino MEC interaction initial "
		       << initial_nucleon_cluster->Pdg() << " final " << final_nucleon_cluster_pdg;
    // Replace this crude fail with the right one for GENIE.
    assert(false);
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

    // remember we pulled this from the prevailing Fermi/SF model, not hard-coded.
    //double ebind1 = removalenergy1; // -0.0250;
    //double ebind2 = removalenergy2; // -0.0250;
    if(removalenergy1 < 0.0 || removalenergy2 < 0.0)
      LOG("MEC", pERROR) << "Removal energy was negative, fail assert " << removalenergy1 << " " << removalenergy2;
    assert(removalenergy1 >= 0.0 && removalenergy2 >= 0.0);

    // cast the removal energy as the energy component of a 4-vector for later.
    TLorentzVector tLVebind(0., 0., 0., -1.0 * (removalenergy1 + removalenergy2));
    
    // RIK you might ask why is this the right place to subtract ebind ?
    // its okay.  physically, I'm subtracting it from q.
    // the energy transfer to the nucleon is 50 MeV less.
    // the energy cost to put this cluster on-shell.

    // what Jan says he does in PRC.86.015504 is this
    // The nucleons are assumed to be in a potential well
    // V = Efermi + 8 MeV.
    // The Fermi energy is subtracted from each initial-state nucleon
    // (I guess he does it dynamically based on Ef = p2/2M or so 
    //     which is what we are doing above, on average)
    // then after the reaction, another 8 MeV is subtracted
    // at that point, a small adjustment to the momentum is needed
    // to keep the nucleon on shell.

    // Did I fix this worry, does this ~50 MeV need to go somewhere in the final state?
    // hmm, but now there is 66 MeV of energy that need to go somewhere.
    // what does GENIE do for QE events in this situation?   

    
    // calculate recoil nucleon cluster 4-momentum
    // tLVebind is intrinsically negative, as it was before.
    p4final_cluster = p4initial_cluster + Q4 + tLVebind;

    
    // Test if the resulting four-vector corresponds to a high-enough invariant mass.
    // Fail the accept if we couldn't put this thing on-shell.
    
    if(p4final_cluster.M() < PDGLibrary::Instance()->Find( final_nucleon_cluster_pdg )->Mass()){
      accept = false;

      // what kinematics causes this
      /*
      std::cout << "RIK Q4 " << iter << " " << Q4.E() << " " << Q4.P() << " xyz " << Q4.Px() << " " << Q4.Py() << " " << Q4.Pz() 
		<< std::endl;

      std::cout << "RIK reject " << iter << " "  << p4final_cluster.Px() << " "
		<< p4final_cluster.Py() << " "
		<< p4final_cluster.Pz() << " "
		<< p4final_cluster.E() << " "
		<< p4final_cluster.M() << " "
		<< PDGLibrary::Instance()->Find( final_nucleon_cluster_pdg )->Mass() << " "
		<< final_nucleon_cluster_pdg 
		<< std::endl;
      */
      
    } else {
      accept = true;
    }

  }  // end accept loop

  // we got here if we accepted the final state two-nucleon system
  // so now we need to write everything to ghep

  // First the initial state nucleons.
  initial_nucleon_cluster->SetMomentum(p4initial_cluster);

  //event->Summary()->InitStatePtr()->TgtPtr()->SetHitNucP4(p4initial_cluster);


  // and the remnant nucleus
  double Mi  = PDGLibrary::Instance()->Find(target_nucleus->Pdg() )-> Mass(); // initial nucleus mass
  remnant_nucleus->SetMomentum(-1.0*p4initial_cluster.Px(),
			       -1.0*p4initial_cluster.Py(),
			       -1.0*p4initial_cluster.Pz(),
			       Mi - p4initial_cluster.E());


  // Now the final nucleon cluster.
  
  // Getting this v4 is a position, i.e. a position within the nucleus (?)
  // possibly it takes on a non-trivial value only for local Fermi gas
  // or for sophisticated treatments of intranuclear rescattering.
  TLorentzVector v4(*neutrino->X4());
  
  // Now write the (undecayed) final two-nucleon system
  GHepParticle p1(final_nucleon_cluster_pdg, kIStDecayedState, 2, -1, -1, -1, p4final_cluster, v4);
  p1.SetRemovalEnergy(removalenergy1 + removalenergy2);
  // actually, this is not an status1 particle, so it is not picked up by the aggregator.
  // and anyway, the aggregator does not run until the very end.
  event->AddParticle(p1);

  interaction->KinePtr()->SetHadSystP4(p4final_cluster);

  LOG("MEC",pDEBUG) << "Generate Initial Hadrons - Done";

}

/*
//_________________________________________________________________________
void MECTensorGenerator::RecoilNucleonCluster    (GHepRecord * event) const
{
    // RIK removed this old method.
}
*/


//_________________________________________________________________________
void MECTensorGenerator::DecayNucleonCluster  (GHepRecord * event) const
{
  LOG("MEC",pDEBUG) << "Decay Nucleon Cluster - Start";

  // Perform a phase-space decay of the nucleon cluster and add its decay
  // products in the event record

  Interaction * interaction = event->Summary();

  // get di-nucleon cluster
  int nucleon_cluster_id = 5;
  GHepParticle * nucleon_cluster = event->Particle(nucleon_cluster_id);
  assert(nucleon_cluster);

  // get decay products
  PDGCodeList pdgv = this->NucleonClusterConstituents(nucleon_cluster->Pdg());
  LOG("MEC", pINFO) << "Decay product IDs: " << pdgv;

  // Get the decay product masses
  vector<int>::const_iterator pdg_iter;
  int iii = 0;
  double * mass = new double[pdgv.size()];
  double   sum  = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
    int pdgc = *pdg_iter;
    double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
    mass[iii++] = m;
    sum += m;
  }

  TLorentzVector * p4d = nucleon_cluster->GetP4();
  TLorentzVector * v4d = nucleon_cluster->GetX4();

  // Rik testing, comment this out or delete it when I'm done
  // this is for checking the action of fPhaseSpaceGenerator 
  // p4d->SetPx(0.0); p4d->SetPy(0.0); p4d->SetPz(0.0);


  // Set the decay
  bool permitted = fPhaseSpaceGenerator.SetDecay(*p4d, pdgv.size(), mass);
  if(!permitted) {
    LOG("MEC",pWARN) << "RIK fPhaseSpaceGenerator says this is not permitted.";
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

  RandomGen * rnd = RandomGen::Instance();
  bool accept_decay=false;
  unsigned int itry=0;
  while(!accept_decay) 
  {
     itry++;

     if(itry > controls::kMaxUnweightDecayIterations) {
       LOG("MEC",pWARN) << "RIK fPhaseSpaceGenerator couldn't do this either.";
       // clean up
       delete [] mass;
       delete p4d;
       delete v4d;
       // throw exception
       LOG("MEC", pWARN)
	 << "Couldn't decay after " 
	 << itry << " iterations";
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select decay after N attempts");
       exception.SwitchOnFastForward();
       throw exception;
     }
     double w  = fPhaseSpaceGenerator.Generate();   
     if(w > wmax) {
       //RIK Jackie deleted a warning here...
       LOG("MEC", pWARN) << " w > wmax " ;
     }
     double gw = wmax * rnd->RndDec().Rndm();
     accept_decay = (gw<=w);

  } //!accept_decay

  // Insert the decay products in the event record
  // The default kPhaseSpaceGenerator class, TGenPhaseSpace,
  // does the boost to lab frame itself.
  TLorentzVector v4(*v4d); 
  GHepStatus_t ist = kIStHadronInTheNucleus;
  int idp = 0;
  double sumE = 0.0;
  double sumpx = 0.0;  double sumpy = 0.0;  double sumpz = 0.0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
     int pdgc = *pdg_iter;
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(idp);
     event->AddParticle(pdgc, ist, nucleon_cluster_id,-1,-1,-1, *p4fin, v4);
     idp++;
     sumE += p4fin->E();
     sumpx += p4fin->Px();
     sumpy += p4fin->Py();
     sumpz += p4fin->Pz();
  }
  
  // do I SetW here?  Do I need to ?
  double W2 = sumE*sumE - sumpx*sumpx - sumpy*sumpy - sumpz*sumpz;
  if(W2 < 0.0){
    LOG("MEC", pERROR) << "Negative W2 " << W2;
    W2 = 0.0;
  }
  interaction->KinePtr()->SetW(TMath::Sqrt(W2), true);

  event->AddParticle(kPdgBindino, kIStStableFinalState,
		     -1,-1,-1,-1, 0,0,0,nucleon_cluster->RemovalEnergy(), 0,0,0,0);

  // Clean-up
  delete [] mass;
  delete p4d;
  delete v4d;

  LOG("MEC",pDEBUG) << "Decay Nucleon Cluster - Done";
}
//_________________________________________________________________________
PDGCodeList MECTensorGenerator::NucleonClusterConstituents(int pdgc) const
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
void MECTensorGenerator::Configure(const Registry & config)   
{
  Algorithm::Configure(config);
  this->LoadConfig();
} 
//___________________________________________________________________________ 
void MECTensorGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void MECTensorGenerator::LoadConfig(void)
{
  fNuclModel = 0;
      
  RgKey nuclkey = "NuclearModel";
  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);
}
//___________________________________________________________________________


