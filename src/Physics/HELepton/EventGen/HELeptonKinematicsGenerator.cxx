//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC & Harvard University
*/
//____________________________________________________________________________

#include "Physics/HELepton/EventGen/HELeptonKinematicsGenerator.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
HELeptonKinematicsGenerator::HELeptonKinematicsGenerator() :
KineGeneratorWithCache("genie::HELeptonKinematicsGenerator")
{

}
//___________________________________________________________________________
HELeptonKinematicsGenerator::HELeptonKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::HELeptonKinematicsGenerator", config)
{

}
//___________________________________________________________________________
HELeptonKinematicsGenerator::~HELeptonKinematicsGenerator()
{

}
//___________________________________________________________________________
void HELeptonKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("HELeptonKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  Interaction * interaction = evrec->Summary();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = this->MaxXSec(evrec);

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(proc_info.IsPhotonCoherent()) {

    double nupdg = interaction->InitState().ProbePdg(); 

    double n2min =  0.;
    double n2max =  1.;
    double n3min =  0.;
    double n3max =  1.;
    double dn2   = n2max-n2min;
    double dn3   = n3max-n3min;

    double n1max = 0.;
    double n1min = 0.;
    if      (pdg::IsNuE  (TMath::Abs(nupdg))) { n1min = 1.-fDeltaN1NuE;   n1max = 1.+fDeltaN1NuE; }
    else if (pdg::IsNuMu (TMath::Abs(nupdg))) { n1min = 1.-fDeltaN1NuMu;  n1max = 1.+fDeltaN1NuMu; }
    else if (pdg::IsNuTau(TMath::Abs(nupdg))) { n1min = 1.-fDeltaN1NuTau; n1max = 1.+fDeltaN1NuTau; }
    double dn1 = n1max-n1min;


    //-- Try to select a valid inelastisity y
    double xsec = -1;
    unsigned int iter = 0;
    bool accept = false;
    
    while(1) {
      iter++;
      if(iter > 1000000) {
        LOG("HELeptonKinematics", pWARN)
              << "*** Could not select a valid y after "
                                              << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
      }

      double n2 = n2min + dn2 * rnd->RndKine().Rndm();
      double n3 = n3min + dn3 * rnd->RndKine().Rndm();
      double n1 = n1min + dn1 * rnd->RndKine().Rndm();
      n1 = (n1>1.) ? n1-2. : n1;

      interaction->KinePtr()->SetKV(kKVn1,n1);
      interaction->KinePtr()->SetKV(kKVn2,n2);
      interaction->KinePtr()->SetKV(kKVn3,n3);

      LOG("HELeptonKinematics", pDEBUG) << "Trying: n1 = " << n1 << ", n2 = " << n2 << ", n3 = " << n3;

      //-- computing cross section for the current kinematics
      xsec = fXSecModel->XSec(interaction, kPSn1n2n3fE);

      this->AssertXSecLimits(interaction, xsec, xsec_max);

      double t = xsec_max * rnd->RndKine().Rndm();
      LOG("HELeptonKinematics", pDEBUG) << "xsec= "<< xsec<< ", J= 1, Rnd= "<< t;

      accept = (t<xsec);

      //-- If the generated kinematics are accepted, finish-up module's job
      if(accept) {
        LOG("HELeptonKinematics", pINFO) << "Selected: n1 = " << n1 << ", n2 = " << n2 << ", n3 = " << n3;

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec,kPSn1n2n3fE);

        // lock selected kinematics & clear running values
        interaction->KinePtr()->ClearRunningValues();

        return;
      }
    }// iterations

  }
  else {
    double n1min = -1.;
    double n1max =  1.;
    double n2min =  0.;
    double n2max =  1.;
    double dn1   = n1max-n1min;
    double dn2   = n2max-n2min;

    //-- Try to select a valid inelastisity y
    double xsec = -1;
    unsigned int iter = 0;
    bool accept = false;
    
    while(1) {
      iter++;
      if(iter > 1000000) {
        LOG("HELeptonKinematics", pWARN)
              << "*** Could not select a valid y after "
                                              << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
      }

      double n1 = n1min + dn1 * rnd->RndKine().Rndm();
      double n2 = n2min + dn2 * rnd->RndKine().Rndm();
      interaction->KinePtr()->SetKV(kKVn1,n1);
      interaction->KinePtr()->SetKV(kKVn2,n2);

      LOG("HELeptonKinematics", pDEBUG) << "Trying: n1 = " << n1 << ", n2 = " << n2;

      //-- computing cross section for the current kinematics
      xsec = fXSecModel->XSec(interaction, kPSn1n2fE);

      this->AssertXSecLimits(interaction, xsec, xsec_max);

      double t = xsec_max * rnd->RndKine().Rndm();
      LOG("HELeptonKinematics", pDEBUG) << "xsec= "<< xsec<< ", J= 1, Rnd= "<< t;

      accept = (t<xsec);

      //-- If the generated kinematics are accepted, finish-up module's job
      if(accept) {
        LOG("HELeptonKinematics", pINFO) << "Selected: n1 = " << n1 << ", n2 = " << n2;

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec,kPSn1n2fE);

        // lock selected kinematics & clear running values
        interaction->KinePtr()->ClearRunningValues();

        return;
      }
    }// iterations

  }
}
//___________________________________________________________________________
double HELeptonKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small y step.

  double max_xsec = -1.0;
  
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit");

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(proc_info.IsPhotonCoherent()) {

    utils::gsl::d2Xsec_dn1dn2dn3_E f(fXSecModel,interaction,-1);
    min->SetFunction( f );
    min->SetMaxFunctionCalls(10000);
    min->SetTolerance(0.05);

    min->SetFixedVariable   ( 0, "n1",   1.);
    min->SetLimitedVariable ( 1, "n2",   0.,   0.01,  0., 1.);
    min->SetLimitedVariable ( 2, "n3",   0.,   0.01,  0., 1.);
    min->Minimize();
    min->SetFixedVariable   ( 0, "n1",   1.);
    min->SetLimitedVariable ( 1, "n2",   min->X()[1],   0.01,  TMath::Max(0.,min->X()[1]-0.1), TMath::Min(1.,min->X()[1]+0.1));
    min->SetLimitedVariable ( 2, "n3",   min->X()[2],   0.01,  TMath::Max(0.,min->X()[2]-0.1), TMath::Min(1.,min->X()[2]+0.1));
    min->Minimize();
    interaction->KinePtr()->SetKV(kKVn1,min->X()[0]);
    interaction->KinePtr()->SetKV(kKVn2,min->X()[1]);
    interaction->KinePtr()->SetKV(kKVn3,min->X()[2]);
    SLOG("HELeptonKinematics", pDEBUG) << "Minimum found -> n1: 1, n2: " << min->X()[1] << ", n3: " << min->X()[2] << ", xsec: " << fXSecModel->XSec(interaction, kPSn1n2n3fE);
    max_xsec = TMath::Max(fXSecModel->XSec(interaction, kPSn1n2n3fE),max_xsec);

    min->SetFixedVariable   ( 0, "n1",   -1.);
    min->SetLimitedVariable ( 1, "n2",   0.,   0.01,  0., 1.);
    min->SetLimitedVariable ( 2, "n3",   0.,   0.01,  0., 1.);
    min->Minimize();
    min->SetFixedVariable   ( 0, "n1",   -1.);
    min->SetLimitedVariable ( 1, "n2",   min->X()[1],   0.01,  TMath::Max(0.,min->X()[1]-0.1), TMath::Min(1.,min->X()[1]+0.1));
    min->SetLimitedVariable ( 2, "n3",   min->X()[2],   0.01,  TMath::Max(0.,min->X()[2]-0.1), TMath::Min(1.,min->X()[2]+0.1));
    interaction->KinePtr()->SetKV(kKVn1,min->X()[0]);
    interaction->KinePtr()->SetKV(kKVn2,min->X()[1]);
    interaction->KinePtr()->SetKV(kKVn3,min->X()[2]);
    SLOG("HELeptonKinematics", pDEBUG) << "Minimum found -> n1: -1, n2: " << min->X()[1] << ", n3: " << min->X()[2] << ", xsec: " << fXSecModel->XSec(interaction, kPSn1n2n3fE);
    max_xsec = TMath::Max(fXSecModel->XSec(interaction, kPSn1n2n3fE),max_xsec);

  } 
  else {

    const int Nscan = 100;
    const int n1min = -1.;
    const int n1max =  1.;
    const int n2min =  0.;
    const int n2max =  1.;
    const double dn1 = (n1max-n1min)/(double)Nscan;
    const double dn2 = (n2max-n2min)/(double)
    Nscan;

    double scan_n1 = 0.;
    double scan_n2 = 0.;
    for (int i=0; i<Nscan; i++) {
      double n1 = n1min + dn1*i;
      for (int j=0; j<Nscan; j++) {
        double n2 = n2min + dn2*j;
        interaction->KinePtr()->SetKV(kKVn1,n1);
        interaction->KinePtr()->SetKV(kKVn2,n2);
        double dxsec = fXSecModel->XSec(interaction, kPSn1n2fE);
        if ( dxsec > max_xsec ) {
          scan_n1 = n1;
          scan_n2 = n2;
          max_xsec = dxsec;
        }    
      }    
    }

    utils::gsl::d2Xsec_dn1dn2_E f(fXSecModel,interaction,-1);
    min->SetFunction( f );
    min->SetMaxFunctionCalls(10000);
    min->SetTolerance(0.05);
    min->SetLimitedVariable ( 0, "n1", scan_n1, 0.001, TMath::Max(-1.,scan_n1-0.1), TMath::Min(1.,scan_n1+0.1));
    min->SetLimitedVariable ( 1, "n2", scan_n2,   0.1, TMath::Max(-0.,scan_n2-0.1), TMath::Min(1.,scan_n2+0.1));
    min->Minimize();
    interaction->KinePtr()->SetKV(kKVn1,min->X()[0]);
    interaction->KinePtr()->SetKV(kKVn2,min->X()[1]);
    max_xsec = fXSecModel->XSec(interaction, kPSn1n2fE);
    SLOG("HELeptonKinematics", pDEBUG) << "Minimum found -> n1: " << min->X()[0] << ", n2: " << min->X()[1];

  }

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

  SLOG("HELeptonKinematics", pDEBUG) << interaction->AsString();
  SLOG("HELeptonKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("HELeptonKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________
double HELeptonKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void HELeptonKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HELeptonKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HELeptonKinematicsGenerator::LoadConfig(void)
{
// Reads its configuration data from its configuration Registry and loads them
// in private data members to avoid looking up at the Registry all the time.

  //-- Safety factor for the maximum differential cross section
  GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor,  2. ) ;

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
    assert(fMaxXSecDiffTolerance>=0);

  //-- Sampling range for n1 variable
  GetParamDef( "Delta-N1-NuE",   fDeltaN1NuE,   0. ) ;
  GetParamDef( "Delta-N1-NuMu",  fDeltaN1NuMu,  0. ) ;
  GetParamDef( "Delta-N1-NuTau", fDeltaN1NuTau, 0. ) ;


}
//____________________________________________________________________________
