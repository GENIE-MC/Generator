//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - Feb 15, 2009

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 15, 2009 - CA
   This class was first added in version 2.5.1.
 @ Feb 06, 2013 - CA
   When the value of the differential cross-section for the selected kinematics
   is set to the event, set the corresponding KinePhaseSpace_t value too.

*/
//____________________________________________________________________________

#include <cfloat>

#include <TMath.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/GBuild.h"
#include "Conventions/Controls.h"
#include "Conventions/Constants.h"
#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "Diffractive/DFRKinematicsGenerator.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/RunningThreadInfo.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
using namespace genie::utils;

//___________________________________________________________________________
DFRKinematicsGenerator::DFRKinematicsGenerator() :
KineGeneratorWithCache("genie::DFRKinematicsGenerator")
{

}
//___________________________________________________________________________
DFRKinematicsGenerator::DFRKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::DFRKinematicsGenerator", config)
{

}
//___________________________________________________________________________
DFRKinematicsGenerator::~DFRKinematicsGenerator()
{

}
//___________________________________________________________________________
void DFRKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("DFRKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- Get the interaction 
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);

  //-- Get neutrino energy and hit 'nucleon mass' 
  const InitialState & init_state = interaction->InitState();
  double Ev  = init_state.ProbeE(kRfHitNucRest);
  double M   = init_state.Tgt().HitNucP4().M(); // can be off m-shell

  //-- Get the physical W range 
  const KPhaseSpace & kps = interaction->PhaseSpace();

  Range1D_t xl = kps.Limits(kKVx);
  Range1D_t yl = kps.Limits(kKVy);

  LOG("DFRKinematics", pNOTICE) << "x: [" << xl.min << ", " << xl.max << "]";
  LOG("DFRKinematics", pNOTICE) << "y: [" << yl.min << ", " << yl.max << "]";

/*
  Range1D_t W  = kps.Limits(kKVW);
  if(W.max <=0 || W.min>=W.max) {
     LOG("DFRKinematics", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("No available phase space");
     exception.SwitchOnFastForward();
     throw exception;
  }
*/

  assert(xl.min>0 && yl.min>0);

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //-- Try to select a valid (x,y) pair using the rejection method

  double dx = xl.max - xl.min;
  double dy = yl.max - yl.min;
  double gx=-1, gy=-1, gW=-1, gQ2=-1, xsec=-1;

  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
       LOG("DFRKinematics", pWARN)
         << " Couldn't select kinematics after " << iter << " iterations";
       evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select kinematics");
       exception.SwitchOnFastForward();
       throw exception;
     }

     //-- random x,y
     gx = xl.min + dx * rnd->RndKine().Rndm();
     gy = yl.min + dy * rnd->RndKine().Rndm();
     interaction->KinePtr()->Setx(gx);
     interaction->KinePtr()->Sety(gy);
     LOG("DFRKinematics", pNOTICE) 
        << "Trying: x = " << gx << ", y = " << gy;

     //-- compute the cross section for current kinematics
     xsec = fXSecModel->XSec(interaction, kPSxyfE);

     //-- decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);
        double t = xsec_max * rnd->RndKine().Rndm();
	double J = 1;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("DFRKinematics", pDEBUG)
              << "xsec= " << xsec << ", J= " << J << ", Rnd= " << t;
#endif
        accept = (t < J*xsec);
     } 
     else {
        accept = (xsec>0);
     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
         // reset trust bits
         interaction->ResetBit(kISkipProcessChk);
         interaction->ResetBit(kISkipKinematicChk);

         // for uniform kinematics, compute an event weight as
         // wght = (phase space volume)*(differential xsec)/(event total xsec)
         if(fGenerateUniformly) {
            double vol     = kinematics::PhaseSpaceVolume(interaction,kPSxyfE);
            double totxsec = evrec->XSec();
            double wght    = (vol/totxsec)*xsec;
            LOG("DFRKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

            // apply computed weight to the current event weight
            wght *= evrec->Weight();
            LOG("DFRKinematics", pNOTICE) << "Current event wght = " << wght;
            evrec->SetWeight(wght);
         }

         // the DFR cross section should be a triple differential cross section
         // d^2xsec/dxdydt where t is the the square of the 4p transfer to the
         // nucleus. The cross section used for kinematical selection should have
         // the t-dependence integrated out. The t-dependence is of the form
         // ~exp(-bt). Now that the x,y kinematical variables have been selected
         // we can generate a t using the t-dependence as a PDF.
         double Epi   = gy*Ev; // pion energy
         double tmax    = 1.0;
         double tmin    = TMath::Min(tmax, TMath::Power(0.5*kPionMass2/Epi,2.));
         double b     = fBeta;
         double tsum  = (TMath::Exp(-b*tmin) - TMath::Exp(-b*tmax))/b; 
         double rt    = tsum * rnd->RndKine().Rndm();
         double gt    = -1.*TMath::Log(-1.*b*rt + TMath::Exp(-1.*b*tmin))/b;

         LOG("DFRKinematics", pNOTICE)
           << "Selected: t = "<< gt << ", from ["<< tmin << ", "<< tmax << "]";

         // compute W,Q2 for selected x,y
         kinematics::XYtoWQ2(Ev,M,gW,gQ2,gx,gy);

         LOG("DFRKinematics", pNOTICE) 
                << "Selected x,y => W = " << gW << ", Q2 = " << gQ2;

         // set the cross section for the selected kinematics
         evrec->SetDiffXSec(xsec*TMath::Exp(-b*gt), kPSxytfE);

         // lock selected kinematics & clear running values
         interaction->KinePtr()->SetW (gW,  true);
         interaction->KinePtr()->SetQ2(gQ2, true);
         interaction->KinePtr()->Setx (gx,  true);
         interaction->KinePtr()->Sety (gy,  true);
         interaction->KinePtr()->Sett (gt,  true);
         interaction->KinePtr()->ClearRunningValues();
         return;
     }
  } // iterations
}
//___________________________________________________________________________
void DFRKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DFRKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DFRKinematicsGenerator::LoadConfig(void)
{
// Reads its configuration data from its configuration Registry and loads them
// in private data members to avoid looking up at the Registry all the time.

  //-- Safety factor for the maximum differential cross section
  fSafetyFactor = fConfig->GetDoubleDef("MaxXSec-SafetyFactor", 1.25);

  //-- Minimum energy for which max xsec would be cached, forcing explicit
  //   calculation for lower eneries
  fEMin = fConfig->GetDoubleDef("Cache-MinEnergy", 0.8);

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  fMaxXSecDiffTolerance = 
                   fConfig->GetDoubleDef("MaxXSec-DiffTolerance",999999.);
  assert(fMaxXSecDiffTolerance>=0);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  fGenerateUniformly = fConfig->GetBoolDef("UniformOverPhaseSpace", false);

  fBeta = 7;
}
//____________________________________________________________________________
double DFRKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But this needs to be fast - do not use a very fine grid.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DFRKinematics", pDEBUG)
      << "Computing max xsec in allowed phase space";
#endif
  double max_xsec = 0.0;

  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t xl = kps.Limits(kKVx);
  Range1D_t yl = kps.Limits(kKVy);

  int    Ny      = 20;
  int    Nx      = 40;
  double xmin    = xl.min;
  double xmax    = xl.max;
  double ymin    = yl.min;
  double ymax    = yl.max;
  double dx      = (xmax-xmin)/(Nx-1);
  double dy      = (ymax-ymin)/(Ny-1);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DFRKinematics", pDEBUG) 
    << "Searching max. in x [" << xmin << ", " << xmax 
    << "], y [" << ymin << ", " << ymax << "]";
#endif
  double xseclast_y = -1;
  bool increasing_y;

  for(int i=0; i<Ny; i++) {
     double gy = ymin + i*dy;
     interaction->KinePtr()->Sety(gy);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("DFRKinematics", pDEBUG) << "y = " << gy;
#endif
     double xseclast_x = -1;
     bool increasing_x;

     for(int j=0; j<Nx; j++) {
        double gx = xmin + j*dx;
        interaction->KinePtr()->Setx(gx);
        kinematics::UpdateWQ2FromXY(interaction);

        double xsec = fXSecModel->XSec(interaction, kPSxyfE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("DFRKinematics", pINFO) 
                << "xsec(y=" << gy << ", x=" << gx << ") = " << xsec;
#endif
        // update maximum xsec
        max_xsec = TMath::Max(xsec, max_xsec);

        increasing_x = xsec-xseclast_x>=0;
        xseclast_x   = xsec;

/*
        // once the cross section stops increasing, I reduce the step size and
        // step backwards a little bit to handle cases that the max cross section
        // is grossly underestimated (very peaky distribution & large step)
        if(!increasing_x) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
          LOG("DFRKinematics", pDEBUG) 
           << "d2xsec/dxdy|x stopped increasing. Stepping back & exiting x loop";
#endif
          double dxn = dx/(Nxb+1);
          for(int ik=0; ik<Nxb; ik++) {
   	     gx = gx - dxn;
             interaction->KinePtr()->Setx(gx);
             kinematics::UpdateWQ2FromXY(interaction);
             xsec = fXSecModel->XSec(interaction, kPSxyfE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
             LOG("DFRKinematics", pINFO) 
                << "xsec(y=" << gy << ", x=" << gx << ") = " << xsec;
#endif
	  }
          break;
        } // stepping back within last bin
*/

     } // x
     increasing_y = max_xsec-xseclast_y>=0;
     xseclast_y   = max_xsec;
     if(!increasing_y) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("DFRKinematics", pDEBUG) 
           << "d2xsec/dxdy stopped increasing. Exiting y loop";
#endif
       break;
     }
  }// y

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  //  max_xsec *= fSafetyFactor;
  //max_xsec *= ( (Ev<3.0) ? 2.5 : fSafetyFactor);
  max_xsec *= 3;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("DFRKinematics", pDEBUG) << interaction->AsString();
  SLOG("DFRKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("DFRKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;
#endif

  return max_xsec;
}
//___________________________________________________________________________

