//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF (Amsterdam)

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Physics/GlashowResonance/EventGen/GLRESKinematicsGenerator.h"
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

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
GLRESKinematicsGenerator::GLRESKinematicsGenerator() :
KineGeneratorWithCache("genie::GLRESKinematicsGenerator")
{

}
//___________________________________________________________________________
GLRESKinematicsGenerator::GLRESKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::GLRESKinematicsGenerator", config)
{

}
//___________________________________________________________________________
GLRESKinematicsGenerator::~GLRESKinematicsGenerator()
{

}
//___________________________________________________________________________
void GLRESKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("GLRESKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

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
  double xsec_max = this->MaxXSec(evrec);
  
  //-- y range
  const KPhaseSpace & kps = evrec->Summary()->PhaseSpace();
  Range1D_t yl = kps.Limits(kKVy);
  double ymin = yl.min;
  double ymax = yl.max;
  double dy   = ymax-ymin;

  double xsec = -1;
  Interaction * interaction = evrec->Summary();

  //-- Try to select a valid inelastisity y
  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("GLRESKinematics", pWARN)
              << "*** Could not select a valid y after "
                                              << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     double y = ymin + dy * rnd->RndKine().Rndm();
     interaction->KinePtr()->Sety(y);

     LOG("GLRESKinematics", pINFO) << "Trying: y = " << y;

     //-- computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSyfE);

     this->AssertXSecLimits(interaction, xsec, xsec_max);

     double t = xsec_max * rnd->RndKine().Rndm();
     LOG("GLRESKinematics", pDEBUG) << "xsec= "<< xsec<< ", J= 1, Rnd= "<< t;

     accept = (t<xsec);

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("GLRESKinematics", pINFO) << "Selected: y = " << y;

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec,kPSyfE);

        // lock selected kinematics & clear running values
        interaction->KinePtr()->Sety(y, true);
        interaction->KinePtr()->ClearRunningValues();

        return;
     }
  }// iterations
}
//___________________________________________________________________________
double GLRESKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small y step.

  const int N  = 100;

  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t yl = kps.Limits(kKVy);
  const double ymin = yl.min;
  const double ymax = yl.max;

  double max_xsec = -1.0;

  double dy = (ymax-ymin)/(N-1);

  for(int i=0; i<N; i++) {
    double y = ymin + i * dy;
    interaction->KinePtr()->Sety(y);
    double xsec = fXSecModel->XSec(interaction, kPSyfE);

    SLOG("GLRESKinematics", pDEBUG) << "xsec(y = " << y << ") = " << xsec;
    max_xsec = TMath::Max(xsec, max_xsec);

  }//y

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

  SLOG("GLRESKinematics", pDEBUG) << interaction->AsString();
  SLOG("GLRESKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("GLRESKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________
double GLRESKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void GLRESKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESKinematicsGenerator::LoadConfig(void)
{
// Reads its configuration data from its configuration Registry and loads them
// in private data members to avoid looking up at the Registry all the time.

  //-- Safety factor for the maximum differential cross section
  GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor,  2. ) ;

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
    assert(fMaxXSecDiffTolerance>=0);

}
//____________________________________________________________________________
