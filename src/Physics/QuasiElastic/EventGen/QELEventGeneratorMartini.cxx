//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

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
#include "Physics/QuasiElastic/EventGen/QELEventGeneratorMartini.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"

#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/PrintUtils.h"

#include "Framework/EventGen/XSecAlgorithmI.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
using namespace genie::utils;

namespace { // anonymous namespace (file only visibility)
  const double eps = std::numeric_limits<double>::epsilon();
}
//___________________________________________________________________________
QELEventGeneratorMartini::QELEventGeneratorMartini() :
KineGeneratorWithCache("genie::QELEventGeneratorMartini")
{

}
//___________________________________________________________________________
QELEventGeneratorMartini::QELEventGeneratorMartini(string config) :
KineGeneratorWithCache("genie::QELEventGeneratorMartini", config)
{

}
//___________________________________________________________________________
QELEventGeneratorMartini::~QELEventGeneratorMartini()
{

}
//___________________________________________________________________________
void QELEventGeneratorMartini::ProcessEventRecord(GHepRecord * event) const
{
  // If we're working with a free nucleon target, then the SuSAv2 calculation
  // isn't set up to handle it. In these cases, we delegate the work to another
  // EventRecordVisitorI object (likely QELEventGenerator) configured to run as
  // a sub-algorithm.
  if ( !event->Summary()->InitState().Tgt().IsNucleus() ) {
    return fFreeNucleonEventGenerator->ProcessEventRecord( event );
  }

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
void QELEventGeneratorMartini::SelectLeptonKinematics (GHepRecord * event) const
{

  // -- Event Properties -----------------------------//
  Interaction * interaction = event->Summary();
  Kinematics * kinematics = interaction->KinePtr();

  // Choose the appropriate minimum Q^2 value based on the interaction
  // mode (this is important for EM interactions since the differential
  // cross section blows up as Q^2 --> 0)
  double Q2min = genie::controls::kMinQ2Limit; // CC/NC limit
  if ( interaction->ProcInfo().IsEM() ) Q2min = genie::utils::kinematics
    ::electromagnetic::kMinQ2Limit; // EM limit

  // The SuSA 1p1h model kinematics works in a system where
  // the whole nuclear target system has no momentum.
  double Enu = interaction->InitState().ProbeE(kRfLab);

  int NuPDG = interaction->InitState().ProbePdg();
  int TgtPDG = interaction->InitState().TgtPdg();
  // interacton vtx
  TLorentzVector v4(*event->Probe()->X4());
  TLorentzVector tempp4(0.,0.,0.,0.);

  GHepParticle * nucleus = event->TargetNucleus();
  bool have_nucleus = nucleus != 0;

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


  // -- Generate and Test the Kinematics----------------------------------//

  RandomGen * rnd = RandomGen::Instance();
  bool accept = false;
  unsigned int iter = 0;
  unsigned int maxIter = kRjMaxIterations * 1000;

  //e-scat xsecs blow up close to theta=0, MC methods won't work so well...
  // NOTE: SuSAv2 1p1h e-scatting has not been validated yet, use with caution
  if ( NuPDG == 11 ) maxIter *= 1000;

  // Get Max XSec:
  double XSecMax = this->MaxXSec( event );

  LOG("Kinematics", pDEBUG) << "Max XSec = " << XSecMax;

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

      // Calculate what we need to enforce the minimum Q^2 cut
      Q0 = Enu - (T + LepMass);
      Q2 = Q3*Q3 - Q0*Q0;

      // Anti neutrino elastic scattering case - include delta in xsec
      if (!have_nucleus){
        genie::utils::mec::Getq0q3FromTlCostl(T, Costh, Enu, LepMass, Q0, Q3);
        Q3 = sqrt(Q0*Q0+2*kNucleonMass*Q0);
        genie::utils::mec::GetTlCostlFromq0q3(Q0, Q3, Enu, LepMass, T, Costh);
        LOG("QELEvent", pDEBUG) << " Anu elastic case. T, Costh:  " << T << ", " << Costh ;
        LOG("QELEvent", pDEBUG) << " Anu elastic case. Q0, Q3:  " << Q0 << ", " << Q3 ;
      }

      // Don't bother doing hard work if the selected Q3 is greater than Q3Max
      // or if Q^2 is less than the minimum allowed value
      if ( Q3 < fQ3Max && Q2 >= Q2min ){

          kinematics->SetKV(kKVTl, T);
          kinematics->SetKV(kKVctl, Costh);
          LOG("QELEvent", pDEBUG) << " T, Costh, Q2: " << T << ", " << Costh << ", " << Q2;

          double XSec = fXSecModel->XSec(interaction, kPSTlctl);

          // Some debugging if things go wrong here...
          if (XSec > XSecMax) {
              LOG("QELEvent", pDEBUG) << "XSec in cm2 is  " << XSec/(units::cm2);
              LOG("QELEvent", pDEBUG) << "XSec in cm2 /neutron is  " << XSec/(units::cm2*pdg::IonPdgCodeToZ(TgtPDG));
              LOG("QELEvent", pDEBUG) << "XSecMax in cm2 /neutron is  " << XSecMax/(units::cm2*pdg::IonPdgCodeToZ(TgtPDG));
              LOG("QELEvent", pERROR) << " T, Costh, Q2: " << T << ", " << Costh << ", " << Q2;
              LOG("QELEvent", pERROR) << "XSec is > XSecMax for nucleus " << TgtPDG << " "
                                 << XSec << " > " << XSecMax
                                 << " don't let this happen.";
          }
          // decide whether to accept or reject these kinematics
          this->AssertXSecLimits( interaction, XSec, XSecMax );
          accept = XSec > XSecMax*rnd->RndKine().Rndm();
          LOG("QELEvent", pINFO) << "Xsec, Max, Accept: " << XSec << ", "
              << XSecMax << ", " << accept;
              LOG("QELEvent", pDEBUG) << "XSec in cm2 /neutron is  " << XSec/(units::cm2*pdg::IonPdgCodeToZ(TgtPDG));
              LOG("QELEvent", pDEBUG) << "XSecMax in cm2 /neutron is  " << XSecMax/(units::cm2*pdg::IonPdgCodeToZ(TgtPDG));


      }// end if passes q3 test
  } // end main accept-reject loop

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
void QELEventGeneratorMartini::AddTargetNucleusRemnant(GHepRecord * event) const
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
void QELEventGeneratorMartini::GenerateNucleon(GHepRecord * event) const
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
    //Q4.Print();

    // get the target nucleon and nucleus.
    // the remnant nucleus is apparently set, except for its momentum.
    GHepParticle * target_nucleus = event->TargetNucleus();
    bool have_nucleus = (target_nucleus != 0);
    //assert(target_nucleus);
    GHepParticle * initial_nucleon = event->HitNucleon();
    assert(initial_nucleon);
    GHepParticle * remnant_nucleus = event->RemnantNucleus();
    //assert(remnant_nucleus);

    // instantiate an empty local target nucleus, so I can use existing methods
    // to get a momentum from the prevailing Fermi-motion distribution
    int tgtpdg;
    if(have_nucleus) tgtpdg = target_nucleus->Pdg();
    else tgtpdg = kPdgTgtFreeP;
    Target tgt(tgtpdg);

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
    // Possible to produce kinematically unallowed (Pauli blocked) nucleon.
    // Find out, and if so choose them again with this accept/reject loop.
    // Pauli blocking was included in the SuSAv2 tensor tables, so it should not
    // be allowed to affect the inclusive xsec.
    while(!accept){
        iter++;
        if(iter > kRjMaxIterations) {
            // error if try too many times
            LOG("QELEvent", pWARN)
                << "Couldn't select a valid nucleon after "
                << iter << " iterations";
            event->EventFlags()->SetBitNumber(kKineGenErr, true);
            genie::exceptions::EVGThreadException exception;
            exception.SetReason("Couldn't select initial hadron kinematics");
            exception.SwitchOnFastForward();
            throw exception;
        }

        // generate nucleons
        // tgt is a Target object for local use, just waiting to be filled.
        // this sets the pdg of each nucleon and its momentum from user chosen nuclear model

        double hitNucPos = initial_nucleon->X4()->Vect().Mag();
        tgt.SetHitNucPdg(initial_nucleon_pdg);
        fNuclModel->GenerateNucleon(tgt,hitNucPos);
        p3i = fNuclModel->Momentum3();

        // Default: Calculate the removal energy as in Guille's thesis - this is a simplicification of
        // a fairly complex aproach employed in SuSAv2, but we expect it to work pretty well.
        // We should write something about this in the implementation technical paper ...
        // IMPORTANT CAVEAT: By default we choose to allow the binding energy to depend on the interaction
        // (as it should), but this means we don't corrolate the chosen Eb with the intial nucleon
        // momentum. Therefore we can sometimes have initial state nucleons with KE>Eb. This isn't
        // great, but fot the moment it's what we have (working on improvments).

        double q3 = Q4.Vect().Mag();

        if(!have_nucleus){
          // For elastic case no Fermi momentum and no Eb
          removalenergy = 0.0;
          p3i.SetXYZ(0.0,0.0,0.0);
        }
        else if(fForceEbFromModel) removalenergy = fNuclModel->RemovalEnergy();
        else if(fForceFixEb) removalenergy = fEbOR;
        else{
          if(q3<0.827){
            removalenergy = -0.017687 + 0.0564*q3;
          }
          else{
            removalenergy = 0.0289558;
          }
          // The above is for Carbon, should shift for different targets
          removalenergy += (fNuclModel->RemovalEnergy()) - fEbC;
          if(removalenergy<0.005) removalenergy=0.005;
        }

        //removalenergy=0.0;
        //initial_nucleon->SetRemovalEnergy(removalenergy);

        LOG("QELEvent",pDEBUG) << "q3 for this event:" << q3;
        LOG("QELEvent",pDEBUG) << "Removal energy:" << removalenergy;

        // Now write down the initial nucleon four-vector for this choice
        double mass = PDGLibrary::Instance()->Find( initial_nucleon_pdg )->Mass();
        double mass2 = mass*mass;
        double energy = TMath::Sqrt(p3i.Mag2() + mass2);
        p4initial_nucleon.SetPxPyPzE(p3i.Px(),p3i.Py(),p3i.Pz(),energy);

        //One rather unsubtle option for making sure nucleons remain bound.
        //This will give us bound nucleons with a sensible missing eneergy
        //but the distribution of Fermi motion will look crazy.
        // Anything we do is wrong (without semi-inclusive inputs) you just
        // have to decide what is less wrong!
        if(fForceBound && (energy-mass>removalenergy)) continue;

        // cast the removal energy as the energy component of a 4-vector for later.
        TLorentzVector tLVebind(0., 0., 0., -1.0 * (removalenergy));

        // Taken from MEC event generator:
        // calculate recoil nucleon cluster 4-momentum (tLVebind is negative)
        p4final_nucleon = p4initial_nucleon + Q4 + tLVebind;

        // Put on shell as in the Aggregator
        // This is a bit of a horrible approximation but it is hard to think
        // of anything simple and better without semi-inclusive model predictions.
        // However, we are working on an improvments.
        if(have_nucleus){
          double En = p4final_nucleon.E();
          double M = PDGLibrary::Instance()->Find( final_nucleon_pdg )->Mass();
          double pmag_old = p4final_nucleon.Vect().Mag();

          double pmag_new = TMath::Sqrt(utils::math::NonNegative(En*En-M*M));
          double scale    = pmag_new / pmag_old;

          double pxn = scale * p4final_nucleon.Px();
          double pyn = scale * p4final_nucleon.Py();
          double pzn = scale * p4final_nucleon.Pz();

          p4final_nucleon.SetPxPyPzE(pxn,pyn,pzn,En);
        }

        // Extra momentum xfer to keep nucleon on shell - I'll give this to the remnant

        pxb = p4nu.Px()-p4l.Px()-p4final_nucleon.Px();
        pyb = p4nu.Py()-p4l.Py()-p4final_nucleon.Py();
        pzb = p4nu.Pz()-p4l.Pz()-p4final_nucleon.Pz();

        LOG("QELEvent",pDEBUG) << "Remnant momentum is: (" << pxb << ", " << pyb << ", " << pzb << ")";


        // Pauli blocking check:
        // Test if the resulting four-vector corresponds to a high-enough energy.
        // Fail the accept if we couldn't put this thing on-shell.
        // Basically: is energy of the nucleon positive after we subtracting Eb
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

    // we got here if we accepted the final state hadron system
    // so now we need to write everything to ghep

    // First the initial state nucleon.
    initial_nucleon->SetMomentum(p4initial_nucleon);

    // and the remnant nucleus
    if(have_nucleus){
      double Mi  = PDGLibrary::Instance()->Find(target_nucleus->Pdg() )-> Mass();
      remnant_nucleus->SetMomentum(pxb,pyb,pzb, Mi - p4initial_nucleon.E() + removalenergy);
    }

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
void QELEventGeneratorMartini::Configure(const Registry & config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//___________________________________________________________________________
void QELEventGeneratorMartini::Configure(string config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//___________________________________________________________________________
void QELEventGeneratorMartini::LoadConfig(void)
{
    fNuclModel = 0;
    RgKey nuclkey = "NuclearModel";
    fNuclModel = dynamic_cast<const NuclearModelI *> ( this->SubAlg(nuclkey) );
    assert( fNuclModel );

    // Sub-algorithm for generating events on free nucleon targets
    // (not handled by the SuSAv2 calculation)
    fFreeNucleonEventGenerator = dynamic_cast<const EventRecordVisitorI* >(
      this->SubAlg("FreeNucleonEventGenerator") );
    assert( fFreeNucleonEventGenerator );

    //-- Maximum q3 in input hadron tensors
    GetParam( "QEL-Q3Max", fQ3Max ) ;

    //-- Whether to force nucleons to be bound
    GetParam( "QEL-ForceBound", fForceBound) ;

    //-- Whether to force Eb to come from the nuclear model
    GetParam( "QEL-ForceEbFromModel", fForceEbFromModel) ;

    //-- Whether to force some fixed Eb
    GetParam( "QEL-ForceFixedEb", fForceFixEb) ;
    GetParam( "QEL-EbOR", fEbOR) ;

    //-- Carbon Eb - needed for scaling
    this->GetParam( "RFG-NucRemovalE@Pdg=1000060120", fEbC);

    //-- Safety factor for the maximum differential cross section
    GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor , 5.00 ) ;

    //-- Minimum energy for which max xsec would be cached, forcing explicit
    //   calculation for lower energies
    //-- I've set this to an extremely high value to avoid problems with the
    //   SuSAv2 model seen during testing. Everything seems to work OK when
    //   the cache is disabled. - S. Gardiner, 1 July 2020
    GetParamDef( "Cache-MinEnergy", fEMin, 1000.00 ) ;

    //-- Maximum allowed fractional cross section deviation from maxim cross
    //   section used in rejection method
    GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
    assert( fMaxXSecDiffTolerance >= 0. );

    //-- Generate kinematics uniformly over allowed phase space and compute
    //   an event weight? NOT IMPLEMENTED FOR SUSA YET!
    GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false ) ;
}
//____________________________________________________________________________
double QELEventGeneratorMartini::ComputeMaxXSec(
  const Interaction * interaction ) const
{
  // Computes the maximum differential cross section in the requested phase
  // space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
  // method and the value is cached at a circular cache branch for retrieval
  // during subsequent event generation.
  // The computed max differential cross section does not need to be the exact
  // maximum. The number used in the rejection method will be scaled up by a
  // safety factor. But it needs to be fast - do not use very small steps.

  double max_xsec = utils::mec::GetMaxXSecTlctl( *fXSecModel, *interaction );

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  max_xsec *= fSafetyFactor;

  #ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    SLOG("QELKinematics", pDEBUG) << interaction->AsString();
    SLOG("QELKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
    SLOG("QELKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;
  #endif

  return max_xsec;
}
//___________________________________________________________________________
