//____________________________________________________________________________
/*
  Copyright (c) 2003-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

  Author: Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
  University of Sussex

  Costas Andreopoulos <c.andreopoulos \at cern.ch>
  University of Liverpool
*/
//____________________________________________________________________________

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
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/DarkNeutrino/EventGen/COHDNuEventGenerator.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
COHDNuEventGenerator::COHDNuEventGenerator() :
  EventRecordVisitorI("genie::COHDNuEventGenerator")
{

}
//___________________________________________________________________________
COHDNuEventGenerator::COHDNuEventGenerator(string config) :
  EventRecordVisitorI("genie::COHDNuEventGenerator", config)
{

}
//___________________________________________________________________________
COHDNuEventGenerator::~COHDNuEventGenerator()
{

}
//___________________________________________________________________________
void COHDNuEventGenerator::ProcessEventRecord(GHepRecord * event) const
{
  this -> GenerateKinematics        (event);
  this -> AddFinalStateDarkNeutrino (event);
  this -> AddRecoilNucleus          (event);
}
//___________________________________________________________________________
void COHDNuEventGenerator::GenerateKinematics(GHepRecord * event) const
{
  Interaction * interaction = event->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  // Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  // Get the kinematical limits
  const InitialState & init_state = interaction -> InitState();
  double E_nu  = init_state.ProbeE(kRfLab);

  // Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  double gDNuE   = -1; // generated Dark Neutrino Energy
  double gxsec = -1; // dsig/dDENu at generated DNuE

  // Generate kinematics
  if(fGenerateUniformly) {
    LOG("COHDNu", pFATAL)
      << "Option to generate kinematics uniformly not supported";
    exit(1);
  }
  else {
    // For the subsequent kinematic selection with the rejection method:
    // Calculate the max differential cross section.
    // Always at Q^2 = 0 for energies and model tested,
    // but go ahead and do the calculation nevertheless.
    utils::gsl::dXSec_dEDNu_E xsec_func(fXSecModel, interaction, fDNuMass, -1.);
    Range1D_t DNuEnergy = xsec_func.IntegrationRange();
    double dDNuE = DNuEnergy.max - DNuEnergy.min;
    ROOT::Math::BrentMinimizer1D minimizer;
    minimizer.SetFunction( xsec_func, DNuEnergy.min, DNuEnergy.max);
    minimizer.Minimize(1000, 1, 1E-5);
    double DNuE_for_xsec_max = minimizer.XMinimum();
    double xsec_max = -1. * xsec_func(DNuE_for_xsec_max); // xsec in units of 1E-38 cm2/GeV

    // Try to select a valid E_N
    unsigned int iter = 0;
    while(1) {
      iter++;
      if(iter > kRjMaxIterations) {
        LOG("COHDNu", pWARN)
          << "*** Could not select a valid DNuE after " << iter << " iterations";
        event->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
      } // max iterations

      gDNuE = DNuEnergy.min + dDNuE * rnd->RndKine().Rndm();
      LOG("COHDNu", pINFO) << "Trying: E_N = " << gDNuE;

      // Computing cross section for the current kinematics
      gxsec = -1. * xsec_func(gDNuE);

      if(gxsec > xsec_max) {
        double frac = TMath::Abs(gxsec-xsec_max)/xsec_max;
        if(frac > fMaxXSecDiffTolerance) {
          LOG("COHDNu", pWARN)
            << "Current computed cross-section (" << gxsec/(units::cm2)
            << " cm2/GeV^2) exceeds the maximum cross-section ("
            << xsec_max/(units::cm2) << " beyond the specified tolerance";
        }
      }

      // Decide whether to accept the current kinematic point
      double t = fSafetyFactor * xsec_max * rnd->RndKine().Rndm();
      //this->AssertXSecLimits(interaction, gxsec, xsec_max);
      LOG("COHDNu", pINFO)
        << "dxsec/dQ2 = " << gxsec/(units::cm2) << " cm2/GeV^2"
        << "J = 1, rnd = " << t;
      bool accept = (t<gxsec);
      if(accept) break; // exit loop
    } // 1
  } // generate uniformly

  // reset bits
  interaction->ResetBit(kISkipProcessChk);
  interaction->ResetBit(kISkipKinematicChk);

  double M_target = interaction->InitState().Tgt().Mass();

  double ETimesM = E_nu * M_target;
  double EPlusM  = E_nu + M_target;

  double p_DNu = TMath::Sqrt(gDNuE*gDNuE - fDNuMass2);
  double cos_theta_DNu = (gDNuE*EPlusM - ETimesM - 0.5*fDNuMass2) / (E_nu * p_DNu);
  double theta_DNu = TMath::ACos(cos_theta_DNu);

  // Take a unit vector along the neutrino direction @ the LAB
  GHepParticle * probe  = event->Probe();
  TVector3 unit_nudir = probe->P4()->Vect().Unit();

  TVector3 DNu_3vector = TVector3(0,0,0);
  double phi = 2.*kPi * rnd->RndKine().Rndm();
  DNu_3vector.SetMagThetaPhi(p_DNu, theta_DNu, phi);
  DNu_3vector.RotateUz(unit_nudir);
  TLorentzVector P4_DNu = TLorentzVector(DNu_3vector, gDNuE);
  interaction->KinePtr()->SetFSLeptonP4(P4_DNu);

  // lock selected kinematics & clear running values
  double gQ2 = -(P4_DNu - *(probe->P4())).M2();
  interaction->KinePtr()->SetQ2(gQ2, true);
  interaction->KinePtr()->ClearRunningValues();

  // Set the cross section for the selected kinematics
  event->SetDiffXSec(gxsec * (1E-38*units::cm2), kPSEDNufE);
}
//___________________________________________________________________________
void COHDNuEventGenerator::AddFinalStateDarkNeutrino(GHepRecord * event) const
{
  GHepParticle * probe  = event->Probe();

  const TLorentzVector & vtx = *(probe->X4());
  TLorentzVector x4l(vtx);  // position 4-vector

  event->AddParticle(probe -> Pdg() > 0 ? kPdgDarkNeutrino : kPdgAntiDarkNeutrino,
                     kIStDecayedState,
                     event->ProbePosition(), event->TargetNucleusPosition(),
                     -1,-1, event->Summary()->Kine().FSLeptonP4(), x4l);
}
//___________________________________________________________________________
void COHDNuEventGenerator::AddRecoilNucleus(GHepRecord * event) const
{
  GHepParticle * probe  = event->Probe();
  GHepParticle * target = event->TargetNucleus();
  GHepParticle * fsl    = event->Particle(probe->FirstDaughter());

  const TLorentzVector & p4probe  = * ( probe  -> P4() );
  const TLorentzVector & p4target = * ( target -> P4() );
  const TLorentzVector & p4fsl    = * ( fsl    -> P4() );

  const TLorentzVector & p4recoil = p4probe + p4target - p4fsl;

  LOG("COHDNu", pNOTICE)
    << "Recoil 4-momentum: " << utils::print::P4AsString(&p4recoil);

  const TLorentzVector & vtx = *(probe->X4());

  event->AddParticle(event->TargetNucleus()->Pdg(),
                     kIStStableFinalState,
                     event->TargetNucleusPosition(), event->ProbePosition(),
                     -1,-1, p4recoil, vtx);
}
//___________________________________________________________________________
void COHDNuEventGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHDNuEventGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHDNuEventGenerator::LoadConfig(void)
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

  fDNuMass = 0.;
  this->GetParam("Dark-NeutrinoMass", fDNuMass);
  fDNuMass2 = fDNuMass * fDNuMass;

}
//____________________________________________________________________________
