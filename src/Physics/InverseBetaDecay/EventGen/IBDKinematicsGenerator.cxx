//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Corey Reed <cjreed \at nikhef.nl> 
         using code from the QELKinematicGenerator written by
         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/InverseBetaDecay/EventGen/IBDKinematicsGenerator.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
using namespace genie::utils;

//___________________________________________________________________________
IBDKinematicsGenerator::IBDKinematicsGenerator() :
KineGeneratorWithCache("genie::IBDKinematicsGenerator")
{

}
//___________________________________________________________________________
IBDKinematicsGenerator::IBDKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::IBDKinematicsGenerator", config)
{

}
//___________________________________________________________________________
IBDKinematicsGenerator::~IBDKinematicsGenerator()
{

}
//___________________________________________________________________________
void IBDKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("IBD", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- Get the interaction and set the 'trust' bits
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Note: The kinematic generator would be using the free nucleon cross
  //   section (even for nuclear targets) so as not to double-count nuclear
  //   suppression. This assumes that a) the nuclear suppression was turned
  //   on when computing the cross sections for selecting the current event 
  //   and that b) if the event turns out to be unphysical (Pauli-blocked) 
  //   the next attempted event will be forced to QEL again.
  //   (discussion with Hugh - GENIE/NeuGEN integration workshop - 07APR2006
  interaction->SetBit(kIAssumeFreeNucleon);

  //-- Get the limits for the generated Q2
  const KPhaseSpace & kps = interaction->PhaseSpace();
  const Range1D_t Q2 = kps.Limits(kKVQ2);
  
  if(Q2.max <=0 || Q2.min>=Q2.max) {
     LOG("IBD", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("No available phase space");
     exception.SwitchOnFastForward();
     throw exception;
  }

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  const double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //-- Try to select a valid Q2 using the rejection method

  // kinematical limits
  const double Q2min  = Q2.min;
  const double Q2max  = Q2.max;
  double xsec   = -1.;
  double gQ2    =  0.;

  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("IBD", pWARN)
          << "Couldn't select a valid Q^2 after " << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }
     
     //-- Generate a Q2 value within the allowed phase space
     gQ2 = Q2min + (Q2max-Q2min) * rnd->RndKine().Rndm();
     interaction->KinePtr()->SetQ2(gQ2);
     LOG("IBD", pINFO) << "Trying: Q^2 = " << gQ2;

     //-- Computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSQ2fE);

     //-- Decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);
        const double t = xsec_max * rnd->RndKine().Rndm();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("IBD", pDEBUG)
            << "dxsec/dQ2 = " << xsec << ", rnd = " << t;
#endif
        accept = (t < xsec);
     } else {
        accept = (xsec>0);
     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("IBD", pINFO) << "Selected: Q^2 = " << gQ2;

        // reset bits
        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);
        interaction->ResetBit(kIAssumeFreeNucleon);

        // compute the rest of the kinematical variables

        // get neutrino energy at struck nucleon rest frame and the
        // struck nucleon mass (can be off the mass shell)
        const InitialState & init_state = interaction->InitState();
        const double E = init_state.ProbeE(kRfHitNucRest);
        const double M = init_state.Tgt().HitNucP4().M();

        LOG("IBD", pNOTICE) << "E = " << E << ", M = "<< M;

        // The hadronic inv. mass is equal to the recoil nucleon on-shell mass.
        const int rpdgc = interaction->RecoilNucleonPdg();
        assert(rpdgc);
        const double gW = PDGLibrary::Instance()->Find(rpdgc)->Mass();

        LOG("IBD", pNOTICE) << "Selected: W = "<< gW;

        // (W,Q2) -> (x,y)
        double gx=0, gy=0;
        kinematics::WQ2toXY(E,M,gW,gQ2,gx,gy);

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec,kPSQ2fE);

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
        if(fGenerateUniformly) {
	   const double vol     = kinematics::PhaseSpaceVolume(interaction,kPSQ2fE);
	   const double totxsec = evrec->XSec();
	         double wght    = (vol/totxsec)*xsec;
	   LOG("IBD", pNOTICE)  << "Kinematics wght = "<< wght;

	   // apply computed weight to the current event weight
	   wght *= evrec->Weight();
	   LOG("IBD", pNOTICE) << "Current event wght = " << wght;
	   evrec->SetWeight(wght);
        }

        // lock selected kinematics & clear running values
        interaction->KinePtr()->SetQ2(gQ2, true);
        interaction->KinePtr()->SetW (gW,  true);
        interaction->KinePtr()->Setx (gx,  true);
        interaction->KinePtr()->Sety (gy,  true);
        interaction->KinePtr()->ClearRunningValues();

        return;
     }
  }// iterations
}
//___________________________________________________________________________
void IBDKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void IBDKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void IBDKinematicsGenerator::LoadConfig(void)
{
	//-- Safety factor for the maximum differential cross section
	GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor, 1.25 ) ;

	//-- Minimum energy for which max xsec would be cached, forcing explicit
	//   calculation for lower eneries
	GetParamDef( "Cache-MinEnergy", fEMin, 1.00 ) ;

	//-- Maximum allowed fractional cross section deviation from maxim cross
	//   section used in rejection method
	GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
	assert(fMaxXSecDiffTolerance>=0);

	//-- Generate kinematics uniformly over allowed phase space and compute
	//   an event weight?
	GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false ) ;

}
//____________________________________________________________________________
double IBDKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small dQ2 step.

  double max_xsec = 0.0;

  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t rQ2 = kps.Limits(kKVQ2);
  if(rQ2.min <=0 || rQ2.max <= rQ2.min) return 0.;

  //const double logQ2min = TMath::Log(rQ2.min + kASmallNum);
  //const double logQ2max = TMath::Log(rQ2.max - kASmallNum);
  const double logQ2min = TMath::Log(rQ2.min);
  const double logQ2max = TMath::Log(rQ2.max);

  const int N  = 15;
  const int Nb = 10;

  double dlogQ2   = (logQ2max - logQ2min) /(N-1);
  double xseclast = -1;
  bool   increasing;

  for(int i=0; i<N; i++) {
     double Q2 = TMath::Exp(logQ2min + i * dlogQ2);
     interaction->KinePtr()->SetQ2(Q2);
     double xsec = fXSecModel->XSec(interaction, kPSQ2fE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("IBD", pDEBUG)  << "xsec(Q2= " << Q2 << ") = " << xsec;
#endif
     max_xsec = TMath::Max(xsec, max_xsec);
     increasing = xsec-xseclast>=0;
     xseclast   = xsec;

     // once the cross section stops increasing, I reduce the step size and
     // step backwards a little bit to handle cases that the max cross section
     // is grossly underestimated (very peaky distribution & large step)
     if(!increasing) {
       dlogQ2/=(Nb+1);
       for(int ib=0; ib<Nb; ib++) {
	 Q2 = TMath::Exp(TMath::Log(Q2) - dlogQ2);
         if(Q2 < rQ2.min) continue;
         interaction->KinePtr()->SetQ2(Q2);
         xsec = fXSecModel->XSec(interaction, kPSQ2fE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
         LOG("IBD", pDEBUG)  << "xsec(Q2= " << Q2 << ") = " << xsec;
#endif
         max_xsec = TMath::Max(xsec, max_xsec);
       }
       break;
     }
  }//Q^2

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  max_xsec *= fSafetyFactor;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("IBD", pDEBUG) << interaction->AsString();
  SLOG("IBD", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("IBD", pDEBUG) << "Computed using alg = " << *fXSecModel;
#endif

  return max_xsec;
}
//___________________________________________________________________________

