//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - July 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Conventions/Controls.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/IMDKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
IMDKinematicsGenerator::IMDKinematicsGenerator() :
KineGeneratorWithCache("genie::IMDKinematicsGenerator")
{

}
//___________________________________________________________________________
IMDKinematicsGenerator::IMDKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::IMDKinematicsGenerator", config)
{

}
//___________________________________________________________________________
IMDKinematicsGenerator::~IMDKinematicsGenerator()
{

}
//___________________________________________________________________________
void IMDKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("IMDKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //-- y range
  double ymin = kMinY; // the xsec algorithm would internally compute the
  double ymax = kMaxY; // kinematically allowd range, and return 0 if outside
  double dy   = ymax-ymin;

  double xsec = -1;
  Interaction * interaction = evrec->Summary();

  //-- Try to select a valid inelastisity y
  register unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("IMDKinematics", pWARN)
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

     LOG("IMDKinematics", pINFO) << "Trying: y = " << y;

     //-- computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSyfE);

     //-- decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);

        double t = xsec_max * rnd->RndKine().Rndm();
        LOG("IMDKinematics", pDEBUG) << "xsec= "<< xsec<< ", J= 1, Rnd= "<< t;

        accept = (t<xsec);
     } else {
       accept = (xsec>0);
     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("IMDKinematics", pINFO) << "Selected: y = " << y;

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec);

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
         if(fGenerateUniformly) {
           double vol     = kinematics::PhaseSpaceVolume(interaction,kPSyfE);
           double totxsec = evrec->XSec();
           double wght    = (vol/totxsec)*xsec;
           LOG("IMDKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

           // apply computed weight to the current event weight
           wght *= evrec->Weight();
           LOG("IMDKinematics", pNOTICE) << "Current event wght = " << wght;
           evrec->SetWeight(wght);
        }

        // lock selected kinematics & clear running values
        interaction->KinePtr()->Sety(y, true);
        interaction->KinePtr()->ClearRunningValues();

        return;
     }
  }// iterations
}
//___________________________________________________________________________
double IMDKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small y step.

  const int N  = 10;
  const int Nb =  6;

  const double ymin = kMinY;
  const double ymax = kMaxY;

  double max_xsec = -1.0;

  double dy = (ymax-ymin)/(N-1);
  double xseclast = -1;
  bool   increasing;

  for(int i=0; i<N; i++) {
    double y = ymin + i * dy;
    interaction->KinePtr()->Sety(y);
    double xsec = fXSecModel->XSec(interaction, kPSyfE);

    SLOG("IMDKinematics", pDEBUG) << "xsec(y = " << y << ") = " << xsec;
    max_xsec = TMath::Max(xsec, max_xsec);

    increasing = xsec-xseclast>=0;
    xseclast   = xsec;

    // once the cross section stops increasing, I reduce the step size and
    // step backwards a little bit to handle cases that the max cross section
    // is grossly underestimated (very peaky distribution & large step)
    if(!increasing) {
       dy/=(Nb+1);
       for(int ib=0; ib<Nb; ib++) {
	 y = y-dy;
         if(y<ymin) break;
         interaction->KinePtr()->Sety(y);
         xsec = fXSecModel->XSec(interaction, kPSyfE);
         SLOG("IMDKinematics", pDEBUG) << "xsec(y = " << y << ") = " << xsec;
         max_xsec = TMath::Max(xsec, max_xsec);
       }
       break;
    }
  }//y

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

  SLOG("IMDKinematics", pDEBUG) << interaction->AsString();
  SLOG("IMDKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("IMDKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________
double IMDKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void IMDKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void IMDKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void IMDKinematicsGenerator::LoadConfig(void)
{
  fXSecModel = 
      dynamic_cast<const XSecAlgorithmI *> (this->SubAlg("DiffXSecAlg"));
  assert(fXSecModel);

  fSafetyFactor = fConfig->GetDoubleDef("MaxXSec-SafetyFactor", 1.25);
  fEMin         = fConfig->GetDoubleDef("Cache-MinEnergy",     -1.00);

  fMaxXSecDiffTolerance = fConfig->GetDoubleDef("MaxXSec-DiffTolerance",0.);
  assert(fMaxXSecDiffTolerance>=0);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  fGenerateUniformly = fConfig->GetBoolDef("UniformOverPhaseSpace", false);
}
//____________________________________________________________________________
