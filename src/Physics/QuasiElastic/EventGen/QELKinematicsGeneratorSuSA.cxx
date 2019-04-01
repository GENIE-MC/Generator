//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Stephen Dolan

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/QuasiElastic/EventGen/QELKinematicsGeneratorSuSA.h"

#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/PrintUtils.h"

#include "Framework/EventGen/XSecAlgorithmI.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
using namespace genie::utils;

//___________________________________________________________________________
QELKinematicsGeneratorSuSA::QELKinematicsGeneratorSuSA() :
    EventRecordVisitorI("genie::QELKinematicsGeneratorSuSA")
{

}
//___________________________________________________________________________
QELKinematicsGeneratorSuSA::QELKinematicsGeneratorSuSA(string config) :
    EventRecordVisitorI("genie::QELKinematicsGeneratorSuSA", config)
{

}
//___________________________________________________________________________
QELKinematicsGeneratorSuSA::~QELKinematicsGeneratorSuSA()
{

}
//___________________________________________________________________________
void QELKinematicsGeneratorSuSA::ProcessEventRecord(GHepRecord * event) const
{
  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();
  //event->Print();
  this -> SelectLeptonKinematics(event);
  //event->Print();
  this -> AddTargetNucleusRemnant(event);
  //event->Print();
  this -> GenerateNucleon(event);
  //event->Print();
}
//___________________________________________________________________________
void QELKinematicsGeneratorSuSA::SelectLeptonKinematics (GHepRecord * event) const
{

  // -- limit the maximum XS for the accept/reject loop -- //
  // 
  // MaxXSec parameters.  This whole calculation could be in it's own function?
  // these need to lead to a number that is safely large enough, or crash the run.
  double XSecMaxPar0 =  97.8152 ;
  double XSecMaxPar1 =  -2.70978;
  double XSecMaxPar2 =  0.183336;
  double XSecMaxPar3 =  103.455 ;
  double XSecMaxPar4 =  1.51078 ;
  double XSecMaxPar5 =  4.0217  ;
  double XSecOffset  =  6.0;

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
      LOG("QELEvent", pERROR) << "Trying to scale XSecMax for larger nuclei, but "
          << TgtPDG << " isn't a nucleus?";
      assert(false);
    }
  }
  
  // -- Generate and Test the Kinematics----------------------------------//

  RandomGen * rnd = RandomGen::Instance();
  bool accept = false;
  unsigned int iter = 0;
  unsigned int maxIter = kRjMaxIterations*1000;

  //e-scat xsecs blow up close to theta=0, MC methods won't work ...
  if (NuPDG==11){
    maxIter *= 100000;
    CosthMax=0.995;
  }


  // loop over different (randomly) selected T and Costh
  while (!accept) {
      iter++;
      if(iter > maxIter) {
          // error if try too many times
          LOG("QELEvent", pERROR)
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
          // AND set the chosen hadronic system

          /*
          // Taken from Guille's function
          double XSecMax = XSecOffset + XSecMaxPar0 + XSecMaxPar1*Enu + XSecMaxPar2*Enu*Enu
                          + XSecMaxPar3*exp(-XSecMaxPar4*Enu - XSecMaxPar5*Enu*Enu);    
          if(Enu<0.045) XSecMax = 70.;
          //Above is d2sigma/dq0dq3, swap to the dimensions of our xsec calculation: d2sigma/dTmudCthmu
          XSecMax*=(Enu*Plep)/(Q3);
          //Change to per atom and then use GENIE units
          XSecMax*=pdg::IonPdgCodeToZ(TgtPDG);
          XSecMax=XSecMax*units::cm2/1E+39;
          */

          //Old version: 
          double XSecMax = 120.0 * TMath::Power(10.0, 2.2504 * TMath::Log10(Enu) - 9.41158);


          if (NuclearA > 12) XSecMax *=  NuclearAfactorXSecMax;  // Scale it by A, precomputed above.

          //For e-scattering need to change the order of magnitude
          if (NuPDG==11) XSecMax *= 10e10;

          LOG("QELEvent", pDEBUG) << " T, Costh: " << T << ", " << Costh ;

          double XSec = fXSecModel->XSec(interaction, kPSTlctl);

          if (XSec > XSecMax) {
              LOG("QELEvent", pDEBUG) << "XSec in cm2 is  " << XSec/(units::cm2);
              LOG("QELEvent", pDEBUG) << "XSec in cm2 /neutron is  " << XSec/(units::cm2*pdg::IonPdgCodeToZ(TgtPDG));
              LOG("QELEvent", pDEBUG) << "XSecMax in cm2 /neutron is  " << XSecMax/(units::cm2*pdg::IonPdgCodeToZ(TgtPDG));
              LOG("QELEvent", pERROR) << "XSec is > XSecMax for nucleus " << TgtPDG << " " 
                                 << XSec << " > " << XSecMax 
                                 << " don't let this happen.";
          }
          assert(XSec <= XSecMax);
          accept = XSec > XSecMax*rnd->RndKine().Rndm();
          LOG("QELEvent", pINFO) << "Xsec, Max, Accept: " << XSec << ", " 
              << XSecMax << ", " << accept; 
              LOG("QELEvent", pDEBUG) << "XSec in cm2 /neutron is  " << XSec/(units::cm2*pdg::IonPdgCodeToZ(TgtPDG));
              LOG("QELEvent", pDEBUG) << "XSecMax in cm2 /neutron is  " << XSecMax/(units::cm2*pdg::IonPdgCodeToZ(TgtPDG));
          /*
          if(accept){
              // Find out if we should use a p or n initial state
              double myrand = rnd->RndKine().Rndm();

              if (NuPDG > 0) {
                  event->AddParticle(kPdgNeutron, kIStNucleonTarget,
                          1, -1, -1, -1, tempp4, v4);
                  interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgNeutron);
              }
              else {
                  event->AddParticle(kPdgProton, kIStNucleonTarget,
                          1, -1, -1, -1, tempp4, v4);
                  interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgProton); 
              }
              
          } // end if accept
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

  // -- Lepton
  event->AddParticle( pdgc, kIStStableFinalState, momidx, -1, -1, -1, p4l, v4);

  LOG("QELEvent",pDEBUG) << "~~~ LEPTON DONE ~~~";
}
//___________________________________________________________________________
void QELKinematicsGeneratorSuSA::AddTargetNucleusRemnant(GHepRecord * event) const
{
    // add the remnant nuclear target at the GHEP record

    LOG("QELEvent", pINFO) << "Adding final state nucleus";

    double Px = 0;
    double Py = 0;
    double Pz = 0;
    double E  = 0;

    GHepParticle * nucleus = event->TargetNucleus();
    bool have_nucleus = nucleus != 0;
    if (!have_nucleus) return;

    int A = nucleus->A();
    int Z = nucleus->Z();

    int fd = nucleus->FirstDaughter();
    int ld = nucleus->LastDaughter();

    for(int id = fd; id <= ld; id++) {

        // compute A,Z for final state nucleus & get its PDG code and its mass
        GHepParticle * particle = event->Particle(id);
        assert(particle);
        int  pdgc = particle->Pdg();
        bool is_p  = pdg::IsProton (pdgc);
        bool is_n  = pdg::IsNeutron(pdgc);

        if (is_p) Z--;
        if (is_p || is_n) A--;

        Px += particle->Px();
        Py += particle->Py();
        Pz += particle->Pz();
        E  += particle->E();

    }//daughters

    TParticlePDG * remn = 0;
    int ipdgc = pdg::IonPdgCode(A, Z);
    remn = PDGLibrary::Instance()->Find(ipdgc);
    if(!remn) {
        LOG("HadronicVtx", pFATAL)
            << "No particle with [A = " << A << ", Z = " << Z
            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
        assert(remn);
    }

    double Mi = nucleus->Mass();  
    Px *= -1;
    Py *= -1;
    Pz *= -1;
    E = Mi-E;

    // Add the nucleus to the event record
    LOG("QELEvent", pINFO)
        << "Adding nucleus [A = " << A << ", Z = " << Z
        << ", pdgc = " << ipdgc << "]";

    int imom = event->TargetNucleusPosition();
    event->AddParticle(
            ipdgc,kIStStableFinalState, imom,-1,-1,-1, Px,Py,Pz,E, 0,0,0,0);

    LOG("QELEvent", pINFO) << "Done";
    LOG("QELEvent", pINFO) << *event;
}
//___________________________________________________________________________
void QELKinematicsGeneratorSuSA::GenerateNucleon(GHepRecord * event) const
{
    // We need a kinematic limits accept/reject loop here, so generating the
    // initial hadrons is combined with generating the recoil hadrons...

    LOG("QELEvent",pDEBUG) << "Generate Nucleon - Start";

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
    Q4.Print();

    // get the target nucleon and nucleus.
    // the remnant nucleus is apparently set, except for its momentum.
    GHepParticle * target_nucleus = event->TargetNucleus();
    assert(target_nucleus);
    GHepParticle * initial_nucleon = event->HitNucleon();
    assert(initial_nucleon);
    GHepParticle * remnant_nucleus = event->RemnantNucleus();
    assert(remnant_nucleus);

    // instantiate an empty local target nucleus, so I can use existing methods
    // to get a momentum from the prevailing Fermi-motion distribution 
    Target tgt(target_nucleus->Pdg());

    // These things need to be saved through to the end of the accept loop.
    bool accept = false;
    TVector3 p3i;
    unsigned int iter = 0;

    int initial_nucleon_pdg = initial_nucleon->Pdg();
    int final_nucleon_pdg = interaction->RecoilNucleonPdg();

    TLorentzVector p4initial_nucleon;
    TLorentzVector p4final_nucleon;
    double removalenergy;

    //remnant kinematic alterations
    double pxb = 0;
    double pyb = 0;
    double pzb = 0;

    //===========================================================================
    // Choose nucleons from the prevailing fermi-motion distribution.
    // Possible to produce kinematically unallowed nucleon.
    // Find out, and if so choose them again with this accept/reject loop.
    // Some kinematics are especially tough (or at least they were for 2p2h...)
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
        double hitNucPos = initial_nucleon->X4()->Vect().Mag();
        tgt.SetHitNucPdg(initial_nucleon_pdg);
        fNuclModel->GenerateNucleon(tgt,hitNucPos);
        p3i = fNuclModel->Momentum3();

        // Calculate the removal energy as in Guille's thesis - this is a simplicification of 
        // a fairly complex aproach employed in SuSAv2, but we expect it to work pretty well. 
        // We should write something about this in the implementation technical paper ... 
        //removalenergy = fNuclModel->RemovalEnergy();
        double q3 = Q4.Vect().Mag();
        if(q3<0.827){
          removalenergy = -0.017687 + 0.0564*q3;
        }
        else{
          removalenergy = 0.0289558;
        }
        if(removalenergy<0.005) removalenergy=0.005;

        //removalenergy=0.0;
        //initial_nucleon->SetRemovalEnergy(removalenergy);

        LOG("QELEvent",pDEBUG) << "q3 for this event:" << q3;
        LOG("QELEvent",pDEBUG) << "Removal energy:" << removalenergy;

        // Now write down the initial nucleon four-vector for this choice
        double mass2 = PDGLibrary::Instance()->Find( initial_nucleon_pdg )->Mass();
        mass2 *= mass2;
        double energy = TMath::Sqrt(p3i.Mag2() + mass2);
        p4initial_nucleon.SetPxPyPzE(p3i.Px(),p3i.Py(),p3i.Pz(),energy);

        // cast the removal energy as the energy component of a 4-vector for later.
        TLorentzVector tLVebind(0., 0., 0., -1.0 * (removalenergy));

        // Taken from MEC event generator:
        // calculate recoil nucleon cluster 4-momentum (tLVebind is negative)
        p4final_nucleon = p4initial_nucleon + Q4 + tLVebind;

        // Put on shell as in the Aggregator
        // This is a bit of a horrible approximation but it is hard to think 
        // of anything better without semi-inclusive model predictions.
        double En = p4final_nucleon.E();
        double M = PDGLibrary::Instance()->Find( final_nucleon_pdg )->Mass();
        double pmag_old = p4final_nucleon.Vect().Mag();

        double pmag_new = TMath::Sqrt(utils::math::NonNegative(En*En-M*M));
        double scale    = pmag_new / pmag_old;

        double pxn = scale * p4final_nucleon.Px();
        double pyn = scale * p4final_nucleon.Py();
        double pzn = scale * p4final_nucleon.Pz();

        p4final_nucleon.SetPxPyPzE(pxn,pyn,pzn,En);

        // Extra momentum xfer to keep nucleon on shell - I'll give this to the remnant

        pxb = p4nu.Px()-p4final_nucleon.Px()-p4l.Px();
        pyb = p4nu.Py()-p4final_nucleon.Py()-p4l.Py();
        pzb = p4nu.Pz()-p4final_nucleon.Pz()-p4l.Pz();

        LOG("QELEvent",pDEBUG) << "Remnant momentum is: (" << pxb << ", " << pyb << ", " << pzb << ")";


        // Test if the resulting four-vector corresponds to a high-enough energy.
        // Fail the accept if we couldn't put this thing on-shell.
        // Basically: is energy of the nucleon positive after we subtracting Eb
        // Note that this is different than the same code in the MECGernerator
        // This checks using the invarient mass of cluster rather than its energy 
        // EDIT: this stage will now always pass as Eb is no longer subtracted here
        if (p4final_nucleon.E() < PDGLibrary::Instance()->Find(final_nucleon_pdg)->Mass()) {
            accept = false;
            LOG("QELEvent",pDEBUG) << "Rejected nucleon, can't be put on-shell";
            LOG("QELEvent",pDEBUG) << "Nucleon invariant mass:" << p4final_nucleon.M();
            LOG("QELEvent",pDEBUG) << "Nucleon real mass:" << PDGLibrary::Instance()->Find(final_nucleon_pdg)->Mass();
            LOG("QELEvent",pDEBUG) << "Nucleon 4 momenutum:";
            //p4final_nucleon.Print();
            LOG("QELEvent",pDEBUG) << "Removal energy:" << removalenergy;
            LOG("QELEvent",pDEBUG) << "Q4 transfer:";
            //Q4.Print();
        }
        else {
            accept = true;
            LOG("QELEvent",pDEBUG) << "Nucleon accepted, Q4 is";
            //Q4.Print();
            LOG("QELEvent",pDEBUG) << "Initial nucleon mass is" << sqrt((p4initial_nucleon.E()*p4initial_nucleon.E())-(p4initial_nucleon.Vect().Mag()*p4initial_nucleon.Vect().Mag()));
            LOG("QELEvent",pDEBUG) << "Final nucleon mass is" << sqrt((p4final_nucleon.E()*p4final_nucleon.E())-(p4final_nucleon.Vect().Mag()*p4final_nucleon.Vect().Mag()));
            LOG("QELEvent",pDEBUG) << "Nucleon invariant mass:" << p4final_nucleon.M();
            LOG("QELEvent",pDEBUG) << "Nucleon real mass:" << PDGLibrary::Instance()->Find(final_nucleon_pdg)->Mass();
            LOG("QELEvent",pDEBUG) << "Nucleon 4 momenutum:";
            //p4final_nucleon.Print();
            LOG("QELEvent",pDEBUG) << "Removal energy:" << removalenergy;
            LOG("QELEvent",pDEBUG) << "Q4 transfer:";
        }

    }  // end accept loop

    // we got here if we accepted the final state two-nucleon system
    // so now we need to write everything to ghep

    // First the initial state nucleons.
    initial_nucleon->SetMomentum(p4initial_nucleon);

    // and the remnant nucleus
    double Mi  = PDGLibrary::Instance()->Find(target_nucleus->Pdg() )-> Mass();
    remnant_nucleus->SetMomentum(pxb,pyb,pzb,
            Mi - p4initial_nucleon.E() + removalenergy);

    // Now the final nucleon.

    // Getting this v4 is a position, i.e. a position within the nucleus (?)
    // possibly it takes on a non-trivial value only for local Fermi gas
    // or for sophisticated treatments of intranuclear rescattering.
    TLorentzVector v4(*neutrino->X4());

    // Now add the final nucleon 

    interaction->KinePtr()->SetHadSystP4(p4final_nucleon);

    GHepStatus_t ist = (tgt.IsNucleus()) ? kIStHadronInTheNucleus : kIStStableFinalState;
    event->AddParticle(interaction->RecoilNucleonPdg(), ist, event->HitNucleonPosition(),-1,-1,-1, interaction->KinePtr()->HadSystP4(), v4);

    LOG("QELEvent",pDEBUG) << "Generate Nucleon - End";

}
//___________________________________________________________________________
void QELKinematicsGeneratorSuSA::Configure(const Registry & config)   
{
    Algorithm::Configure(config);
    this->LoadConfig();
} 
//___________________________________________________________________________ 
void QELKinematicsGeneratorSuSA::Configure(string config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//___________________________________________________________________________
void QELKinematicsGeneratorSuSA::LoadConfig(void)
{
    fNuclModel = 0;
    RgKey nuclkey = "NuclearModel";
    fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
    assert(fNuclModel);

    GetParam( "QEL-Q3Max", fQ3Max ) ;
}
//___________________________________________________________________________
