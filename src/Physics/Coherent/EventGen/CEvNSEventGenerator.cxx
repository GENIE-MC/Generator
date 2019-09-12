//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
         University of Liverpool & STFC Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/GSLMinimizer1D.h>
#include <Math/BrentMinimizer1D.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/Coherent/EventGen/CEvNSEventGenerator.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
CEvNSEventGenerator::CEvNSEventGenerator() :
EventRecordVisitorI("genie::CEvNSEventGenerator")
{

}
//___________________________________________________________________________
CEvNSEventGenerator::CEvNSEventGenerator(string config) :
EventRecordVisitorI("genie::COHElKinematicsGenerator", config)
{

}
//___________________________________________________________________________
CEvNSEventGenerator::~CEvNSEventGenerator()
{

}
//___________________________________________________________________________
void CEvNSEventGenerator::ProcessEventRecord(GHepRecord * event) const
{
  this -> GenerateKinematics    (event);
  this -> AddFinalStateNeutrino (event);
  this -> AddRecoilNucleus      (event);
}
//___________________________________________________________________________
void CEvNSEventGenerator::GenerateKinematics(GHepRecord * event) const
{
  Interaction * interaction = event->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  // Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  // Get the kinematical limits
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t Q2 = kps.Q2Lim();
  assert(Q2.min > 0. && Q2.min < Q2.max);
  const InitialState & init_state = interaction -> InitState();
  double E  = init_state.ProbeE(kRfLab);
  const double Q2min = Q2.min;
  const double Q2max = Q2.max;
  const double dQ2   = Q2max - Q2min;

  // Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  double gQ2   = -1; // generated Q2
  double gxsec = -1; // dsig/dQ2 at generated Q2

  // Generate kinematics
  if(fGenerateUniformly) {
    LOG("CEvNS", pFATAL)
      << "Option to generate kinematics uniformly not supported";
    exit(1);
/*
    gQ2 = Q2min + dQ2 * rnd->RndKine().Rndm();
    LOG("CEvNS", pINFO) << "Trying: Q2 = " << gQ2;
    interaction->KinePtr()->SetQ2(gQ2);

    // Computing cross section for the current kinematics
    gxsec = fXSecModel->XSec(interaction, kPSQ2fE);
    if(gxsec<=0){
      LOG("CEvNS", pWARN)
        << "Non-positive x-section for selected Q2 = " << gQ2 << "GeV^2";
    }

    double weight = 1; // to implement if fGenerateUniformly option is enabled
    LOG("CEvNS", pDEBUG)  << "Kinematics wght = "<< weight;
    // apply computed weight to the current event weight
    weight *= event->Weight();
    LOG("CEvNS", pDEBUG) << "Current event wght = " << weight;
    event->SetWeight(weight);
*/
  }
  else {
    // For the subsequent kinematic selection with the rejection method:
    // Calculate the max differential cross section.
    // Always at Q^2 = 0 for energies and model tested,
    // but go ahead and do the calculation nevertheless.
    ROOT::Math::IBaseFunctionOneDim * xsec_func =
        new utils::gsl::dXSec_dQ2_E(fXSecModel, interaction,-1.);
    ROOT::Math::BrentMinimizer1D minimizer;
    minimizer.SetFunction(*xsec_func,Q2min,Q2max);
    minimizer.Minimize(1000,1,1E-5);
    double Q2_for_xsec_max = minimizer.XMinimum();
    interaction->KinePtr()->SetQ2(Q2_for_xsec_max);
    double xsec_max = fXSecModel->XSec(interaction, kPSQ2fE);
    delete xsec_func;
    LOG("CEvNS", pNOTICE)
      << "Maximizing dsig(Q2;E = " << E << "GeV)/dQ2 gave a value of "
      << xsec_max/(units::cm2) << " cm2/GeV^2 at Q2 = "
      << Q2_for_xsec_max << " GeV^2";

    // Try to select a valid Q2
    unsigned int iter = 0;
    while(1) {
       iter++;
       if(iter > kRjMaxIterations) {
          LOG("CEvNS", pWARN)
            << "*** Could not select a valid Q2 after " << iter << " iterations";
          event->EventFlags()->SetBitNumber(kKineGenErr, true);
          genie::exceptions::EVGThreadException exception;
          exception.SetReason("Couldn't select kinematics");
          exception.SwitchOnFastForward();
          throw exception;
       } // max iterations

       gQ2 = Q2min + dQ2 * rnd->RndKine().Rndm();
       LOG("CEvNS", pINFO) << "Trying: Q2 = " << gQ2;
       interaction->KinePtr()->SetQ2(gQ2);

       // Computing cross section for the current kinematics
       gxsec = fXSecModel->XSec(interaction, kPSQ2fE);

       if(gxsec > xsec_max) {
          double frac = TMath::Abs(gxsec-xsec_max)/xsec_max;
          if(frac > fMaxXSecDiffTolerance) {
            LOG("CEvNS", pWARN)
              << "Current computed cross-section (" << gxsec/(units::cm2)
              << " cm2/GeV^2) exceeds the maximum cross-section ("
              << xsec_max/(units::cm2) << " beyond the specified tolerance";
          }
       }

       // Decide whether to accept the current kinematic point
       double t = fSafetyFactor * xsec_max * rnd->RndKine().Rndm();
       //this->AssertXSecLimits(interaction, gxsec, xsec_max);
       LOG("CEvNS", pINFO)
         << "dxsec/dQ2 = " << gxsec/(units::cm2) << " cm2/GeV^2"
         << "J = 1, rnd = " << t;
       bool accept = (t<gxsec);
       if(accept) break; // exit loop
    } // 1
  } // generate uniformly

  LOG("CEvNS", pNOTICE) << "Selected Q2 = " << gQ2 << " GeV^2";

  // reset bits
  interaction->ResetBit(kISkipProcessChk);
  interaction->ResetBit(kISkipKinematicChk);

  // lock selected kinematics & clear running values
  interaction->KinePtr()->SetQ2(gQ2, true);
  interaction->KinePtr()->ClearRunningValues();

  // Set the cross section for the selected kinematics
  event->SetDiffXSec(gxsec,kPSQ2fE);
}
//___________________________________________________________________________
void CEvNSEventGenerator::AddFinalStateNeutrino(GHepRecord * event) const
{
  GHepParticle * probe  = event->Probe();
  GHepParticle * target = event->TargetNucleus();

  int target_pdgc = target->Pdg();
  double M = PDGLibrary::Instance()->Find(target_pdgc)->Mass(); // units: GeV

  double Ev = probe->E(); // neutrino energy, units: GeV
  double Q2 = event->Summary()->Kine().Q2(true); // selected momentum transfer, units: GeV^2

  const TLorentzVector & vtx = *(probe->X4());
  TLorentzVector x4l(vtx);  // position 4-vector

  // Compute the final state neutino energy

  double y  = Q2/(2*M*Ev);
  double El = (1-y)*Ev;

  LOG("CEvNS", pNOTICE)
    << "Final state neutrino energy: E = " << El << " GeV";

  // Compute the final state neutrino momentum components
  // along and perpendicular to the incoming neutrino direction

  double ml2 = 0;
  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)

  LOG("CEvNS", pNOTICE)
    << "Final state neutrino momentum components: |p//| = "
    << plp << " GeV, [pT] = " << plt << " GeV";

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);

  // Take a unit vector along the neutrino direction @ the LAB
  TVector3 unit_nudir = probe->P4()->Vect().Unit();

  // Rotate lepton momentum vector from the reference frame (x'y'z') where
  // {z':(neutrino direction), z'x':(theta plane)} to the LAB
  TVector3 p3l(pltx,plty,plp);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in the LAB
  TLorentzVector p4l(p3l,El);

  LOG("CEvNS", pNOTICE)
     << "Final state neutrino 4-momentum: " << utils::print::P4AsString(&p4l);

  event->AddParticle(
    probe->Pdg(), kIStStableFinalState, event->ProbePosition(),
    -1,-1,-1, p4l, x4l);

  event->Summary()->KinePtr()->SetFSLeptonP4(p4l);
}
//___________________________________________________________________________
void CEvNSEventGenerator::AddRecoilNucleus(GHepRecord * event) const
{
  GHepParticle * probe  = event->Probe();
  GHepParticle * target = event->TargetNucleus();
  GHepParticle * fsl    = event->Particle(probe->FirstDaughter());

  const TLorentzVector & p4probe  = * ( probe  -> P4() );
  const TLorentzVector & p4target = * ( target -> P4() );
  const TLorentzVector & p4fsl    = * ( fsl    -> P4() );

  const TLorentzVector & p4recoil = p4probe + p4target - p4fsl;

  LOG("CEvNS", pNOTICE)
     << "Recoil 4-momentum: " << utils::print::P4AsString(&p4recoil);

  const TLorentzVector & vtx = *(probe->X4());

  event->AddParticle(
      event->TargetNucleus()->Pdg(),
      kIStStableFinalState,
      event->TargetNucleusPosition(),
      -1,-1,-1, p4recoil, vtx);
}
//___________________________________________________________________________
void CEvNSEventGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CEvNSEventGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CEvNSEventGenerator::LoadConfig(void)
{
  fXSecModel = 0;

  // max xsec safety factor (for rejection method) and min cached energy
  this->GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor, 1.05 ) ;

  // Generate kinematics uniformly over allowed phase space and compute
  // an event weight?
  this->GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false ) ;

  // Maximum allowed fractional cross section deviation from maxim cross
  // section used in rejection method
  this->GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
  assert(fMaxXSecDiffTolerance>=0);
}
//____________________________________________________________________________
