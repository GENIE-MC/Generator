//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the NuE package from its previous location (EVGModules package)
 @ Feb 06, 2013 - CA
   When the value of the differential cross-section for the selected kinematics
   is set to the event, set the corresponding KinePhaseSpace_t value too.

*/
//____________________________________________________________________________

#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Physics/NuElectron/EventGen/NuEKinematicsGenerator.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
NuEKinematicsGenerator::NuEKinematicsGenerator() :
KineGeneratorWithCache("genie::NuEKinematicsGenerator")
{

}
//___________________________________________________________________________
NuEKinematicsGenerator::NuEKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::NuEKinematicsGenerator", config)
{

}
//___________________________________________________________________________
NuEKinematicsGenerator::~NuEKinematicsGenerator()
{

}
//___________________________________________________________________________
void NuEKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("NuEKinematics", pNOTICE)
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
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);
  
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
        LOG("NuEKinematics", pWARN)
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

     LOG("NuEKinematics", pINFO) << "Trying: y = " << y;

     //-- computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSyfE);

     //-- decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);

        double t = xsec_max * rnd->RndKine().Rndm();
        LOG("NuEKinematics", pDEBUG) << "xsec= "<< xsec<< ", J= 1, Rnd= "<< t;

        accept = (t<xsec);
     } else {
       accept = (xsec>0);
     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("NuEKinematics", pINFO) << "Selected: y = " << y;

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec,kPSyfE);

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
         if(fGenerateUniformly) {
           double vol     = kinematics::PhaseSpaceVolume(interaction,kPSyfE);
           double totxsec = evrec->XSec();
           double wght    = (vol/totxsec)*xsec;
           LOG("NuEKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

           // apply computed weight to the current event weight
           wght *= evrec->Weight();
           LOG("NuEKinematics", pNOTICE) << "Current event wght = " << wght;
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
double NuEKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small y step.

  const int N  = 40;
//const int Nb =  6;

  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t yl = kps.Limits(kKVy);
  const double ymin = yl.min;
  const double ymax = yl.max;

  double max_xsec = -1.0;

  double dy = (ymax-ymin)/(N-1);
//double xseclast = -1;
//bool   increasing;

  for(int i=0; i<N; i++) {
    double y = ymin + i * dy;
    interaction->KinePtr()->Sety(y);
    double xsec = fXSecModel->XSec(interaction, kPSyfE);

    SLOG("NuEKinematics", pDEBUG) << "xsec(y = " << y << ") = " << xsec;
    max_xsec = TMath::Max(xsec, max_xsec);
/*
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
         SLOG("NuEKinematics", pDEBUG) << "xsec(y = " << y << ") = " << xsec;
         max_xsec = TMath::Max(xsec, max_xsec);
       }
       break;
    }
*/
  }//y

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

  SLOG("NuEKinematics", pDEBUG) << interaction->AsString();
  SLOG("NuEKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("NuEKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________
double NuEKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void NuEKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuEKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuEKinematicsGenerator::LoadConfig(void)
{
	GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor, 2.00 ) ;
	GetParamDef( "Cache-MinEnergy", fEMin, 1.00 ) ;

	GetParamDef("MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 0. ) ;
	assert(fMaxXSecDiffTolerance>=0);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
	GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false ) ;

}
//____________________________________________________________________________
