//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Physics/Coherent/EventGen/COHElKinematicsGenerator.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
COHElKinematicsGenerator::COHElKinematicsGenerator() :
KineGeneratorWithCache("genie::COHElKinematicsGenerator")
{

}
//___________________________________________________________________________
COHElKinematicsGenerator::COHElKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::COHElKinematicsGenerator", config)
{

}
//___________________________________________________________________________
COHElKinematicsGenerator::~COHElKinematicsGenerator()
{

}
//___________________________________________________________________________
void COHElKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("COHElKinematics", pNOTICE)
       << "Generating kinematics uniformly over the allowed phase space";
  }

  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //-- Get the kinematical limits for the generated x,y
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t y = kps.YLim();
  assert(y.min>0. && y.max>0. && y.min<1. && y.max<1. && y.min<y.max);
  const double ymin = y.min + kASmallNum;
  const double ymax = y.max - kASmallNum;
  const double dy   = ymax - ymin;

  //------ Try to select a valid y
  unsigned int iter = 0;
  bool accept=false;
  double xsec=-1, gy=-1;

  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("COHElKinematics", pWARN)
            << "*** Could not select a valid y after " << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     gy = ymin + dy * rnd->RndKine().Rndm();
     LOG("COHElKinematics", pINFO) << "Trying: y = " << gy;

     interaction->KinePtr()->Sety(gy);

     // computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSyfE);

     //-- decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        double t = xsec_max * rnd->RndKine().Rndm();
        this->AssertXSecLimits(interaction, xsec, xsec_max);
        LOG("COHElKinematics", pINFO) << "xsec= " << xsec << ", J= 1, Rnd= " << t;
        accept = (t<xsec);
     }
     else { 
        accept = (xsec>0);
     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("COHElKinematics", pNOTICE) << "Selected: y = " << gy;

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
        if(fGenerateUniformly) {
          double vol     = y.max-y.min; 
          double totxsec = evrec->XSec();
          double wght    = (vol/totxsec)*xsec;
          LOG("COHElKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

          // apply computed weight to the current event weight
          wght *= evrec->Weight();
          LOG("COHElKinematics", pNOTICE) << "Current event wght = " << wght;
          evrec->SetWeight(wght);
        }

        // reset bits
        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);

        // lock selected kinematics & clear running values
        double Ev = interaction->InitState().ProbeE(kRfLab);
        double x  = 1;
        double Q2 = 2*kNucleonMass*x*gy*Ev;
        interaction->KinePtr()->Setx(x,   true);
        interaction->KinePtr()->Sety(gy,  true);
        interaction->KinePtr()->Sett(0,   true);
        interaction->KinePtr()->SetW(0,   true);
        interaction->KinePtr()->SetQ2(Q2, true);
        interaction->KinePtr()->ClearRunningValues();

        // set the cross section for the selected kinematics
        // TODO: this is kPSxyfE on the dev branch, but I think kPSyfE is right
        evrec->SetDiffXSec(xsec,kPSyfE);

        return;
     }
  }// iterations
}
//___________________________________________________________________________
double COHElKinematicsGenerator::ComputeMaxXSec(const Interaction * in) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.

  SLOG("COHElKinematics", pDEBUG)
          << "Scanning the allowed phase space {K} for the max(dxsec/d{K})";

  double max_xsec = 0.;
  const int N     = 50;

  const KPhaseSpace & kps = in->PhaseSpace();
  Range1D_t yr = kps.YLim();

  const double logymin = TMath::Log10(yr.min);
  const double logymax = TMath::Log10(yr.max);
  const double dlogy   = (logymax - logymin) /(N-1);

  for(int i=0; i<N; i++) {
     double y = TMath::Power(10, logymin+i*dlogy);
     in->KinePtr()->Sety(y);

     double xsec = fXSecModel->XSec(in, kPSyfE);
     LOG("COHElKinematics", pDEBUG)  << "xsec(y= " << y << ") = " << xsec;
     max_xsec = TMath::Max(max_xsec, xsec);
  }//y

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

  SLOG("COHElKinematics", pDEBUG) << in->AsString();
  SLOG("COHElKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("COHElKinematics", pDEBUG) << "Computed using alg = " << fXSecModel->Id();

  return max_xsec;
}
//___________________________________________________________________________
double COHElKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void COHElKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHElKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHElKinematicsGenerator::LoadConfig(void)
{
  //-- max xsec safety factor (for rejection method) and min cached energy
  GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor, 1.6 ) ;
  GetParamDef( "Cache-MinEnergy", fEMin, -1.0 ) ;

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false ) ;

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
  assert(fMaxXSecDiffTolerance>=0);
}
//____________________________________________________________________________

