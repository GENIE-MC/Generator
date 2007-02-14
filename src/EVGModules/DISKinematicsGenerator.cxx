//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cfloat>

#include <TMath.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Controls.h"
#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/RunningThreadInfo.h"
#include "EVGModules/DISKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
DISKinematicsGenerator::DISKinematicsGenerator() :
KineGeneratorWithCache("genie::DISKinematicsGenerator")
{

}
//___________________________________________________________________________
DISKinematicsGenerator::DISKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::DISKinematicsGenerator", config)
{

}
//___________________________________________________________________________
DISKinematicsGenerator::~DISKinematicsGenerator()
{

}
//___________________________________________________________________________
void DISKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("DISKinematics", pNOTICE)
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
  Range1D_t W  = kps.Limits(kKVW);
  if(W.max <=0 || W.min>=W.max) {
     LOG("DISKinematics", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("No available phase space");
     exception.SwitchOnFastForward();
     throw exception;
  }

  Range1D_t xl = kps.Limits(kKVx);
  Range1D_t yl = kps.Limits(kKVy);

  LOG("DISKinematics", pINFO) << "x: [" << xl.min << ", " << xl.max << "]";
  LOG("DISKinematics", pINFO) << "y: [" << yl.min << ", " << yl.max << "]";

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
       LOG("DISKinematics", pWARN)
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
     kinematics::UpdateWQ2FromXY(interaction);

     LOG("DISKinematics", pINFO) 
        << "Trying: x = " << gx << ", y = " << gy 
        << " (W  = " << interaction->KinePtr()->W()  << ","
        << " (Q2 = " << interaction->KinePtr()->Q2() << ")";

     //-- compute the cross section for current kinematics
     xsec = fXSecModel->XSec(interaction, kPSxyfE);

     //-- decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);
        double t = xsec_max * rnd->RndKine().Rndm();
        //double J = kinematics::Jacobian(interaction,kPSxyfE,kPSlogxlogyfE);
	double J = 1;

        LOG("DISKinematics", pDEBUG)
              << "xsec= " << xsec << ", J= " << J << ", Rnd= " << t;
        accept = (t < J*xsec);
     } 
     else {
        accept = (xsec>0);
     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
         LOG("DISKinematics", pNOTICE) 
            << "Selected:  x = " << gx << ", y = " << gy
            << " (W  = " << interaction->KinePtr()->W()  << ","
            << " (Q2 = " << interaction->KinePtr()->Q2() << ")";

         // reset trust bits
         interaction->ResetBit(kISkipProcessChk);
         interaction->ResetBit(kISkipKinematicChk);

         // set the cross section for the selected kinematics
         evrec->SetDiffXSec(xsec);

         // for uniform kinematics, compute an event weight as
         // wght = (phase space volume)*(differential xsec)/(event total xsec)
         if(fGenerateUniformly) {
            double vol     = kinematics::PhaseSpaceVolume(interaction,kPSxyfE);
            double totxsec = evrec->XSec();
            double wght    = (vol/totxsec)*xsec;
            LOG("DISKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

            // apply computed weight to the current event weight
            wght *= evrec->Weight();
            LOG("DISKinematics", pNOTICE) << "Current event wght = " << wght;
            evrec->SetWeight(wght);
         }

         // compute W,Q2 for selected x,y
         kinematics::XYtoWQ2(Ev,M,gW,gQ2,gx,gy);

         LOG("DISKinematics", pNOTICE) 
                        << "Selected x,y => W = " << gW << ", Q2 = " << gQ2;

         // lock selected kinematics & clear running values
         interaction->KinePtr()->SetW (gW,  true);
         interaction->KinePtr()->SetQ2(gQ2, true);
         interaction->KinePtr()->Setx (gx,  true);
         interaction->KinePtr()->Sety (gy,  true);
         interaction->KinePtr()->ClearRunningValues();
         return;
     }
  } // iterations
}
//___________________________________________________________________________
void DISKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISKinematicsGenerator::LoadConfig(void)
{
// Reads its configuration data from its configuration Registry and loads them
// in private data members to avoid looking up at the Registry all the time.

  //-- Safety factor for the maximum differential cross section
  fSafetyFactor = fConfig->GetDoubleDef("MaxXSec-SafetyFactor", 1.25);

  //-- Minimum energy for which max xsec would be cached, forcing explicit
  //   calculation for lower eneries
  fEMin = fConfig->GetDoubleDef("Cache-MinEnergy", 1.0);

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  fMaxXSecDiffTolerance = fConfig->GetDoubleDef("MaxXSec-DiffTolerance",0.);
  assert(fMaxXSecDiffTolerance>=0);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  fGenerateUniformly = fConfig->GetBoolDef("UniformOverPhaseSpace", false);
}
//____________________________________________________________________________
double DISKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But this needs to be fast - do not use a very fine grid.

  LOG("DISKinematics", pDEBUG)<< "Computing max xsec in allowed phase space";

  double max_xsec = 0.0;

  const InitialState & init_state = interaction->InitState();
  double Ev = init_state.ProbeE(kRfHitNucRest);

  double log10Ev = TMath::Log10(Ev);

  const ProcessInfo & proc = interaction->ProcInfo();
  const Target & tgt = init_state.Tgt();

  // guess the approximate position of the maximum differential xsec to avoid
  // an extensive (time-consuming) phase space search

  double xpeak=-1, ypeak=-1, xwindow=-1, ywindow=-1;
  int Ny=100, Nx=100, Nxb=10;

  if(proc.IsWeakCC()) {
    if( pdg::IsProton(tgt.HitNucPdg()) ) {
      xpeak = 0.18 - 0.055 * log10Ev;
      ypeak = 0.95 - 0.670 * log10Ev + 0.145 * TMath::Power(log10Ev,2);
    } else {
      xpeak = 0.26 - 0.060 * log10Ev;
      ypeak = 0.85 - 0.680 * log10Ev + 0.145 * TMath::Power(log10Ev,2);
    }
    xwindow = .2;
    ywindow = .2;
    Ny  = 10;
    Nx  = 10;
    Nxb =  3;
  } else {
    xpeak   = .1;
    ypeak   = .5;
    xwindow = .2;
    ywindow = .48;
    Ny  = 50;
    Nx  = 50;
    Nxb =  5;
  }

  if(tgt.HitQrkIsSet() && tgt.HitSeaQrk()) {
    xpeak    = .1;
    xwindow  = .1;
  }

  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t xl = kps.Limits(kKVx);
  Range1D_t yl = kps.Limits(kKVy);

  double xmin    = TMath::Max(xpeak-xwindow, TMath::Max(xl.min, 5E-3));
  double xmax    = TMath::Min(xpeak+xwindow, xl.max);
  double ymin    = TMath::Max(ypeak-ywindow, TMath::Max(yl.min, 2E-3));
  double ymax    = TMath::Min(ypeak+ywindow, yl.max);
  double logxmin = TMath::Log10(xmin);
  double logxmax = TMath::Log10(xmax);
  double logymin = TMath::Log10(ymin);
  double logymax = TMath::Log10(ymax);
  double dlogx   = (logxmax-logxmin)/(Nx-1);
  double dlogy   = (logymax-logymin)/(Ny-1);

  LOG("DISKinematics", pDEBUG) 
    << "Searching max. in x [" << xmin << ", " << xmax << "], y [" << ymin << ", " << ymax << "]";

  double xseclast_y = -1;
  bool increasing_y;

  for(int i=0; i<Ny; i++) {
    //double gy = ymin + i*dy;
     double gy = TMath::Power(10., logymin + i*dlogy);
     interaction->KinePtr()->Sety(gy);

     LOG("DISKinematics", pDEBUG) << "y = " << gy;

     double xseclast_x = -1;
     bool increasing_x;

     for(int j=0; j<Nx; j++) {
        //double gx = xmin + j*dx;
	double gx = TMath::Power(10., logxmin + j*dlogx);
        interaction->KinePtr()->Setx(gx);
        kinematics::UpdateWQ2FromXY(interaction);

        double xsec = fXSecModel->XSec(interaction, kPSxyfE);
        LOG("DISKinematics", pINFO) 
                << "xsec(y=" << gy << ", x=" << gx << ") = " << xsec;

        // update maximum xsec
        max_xsec = TMath::Max(xsec, max_xsec);

        increasing_x = xsec-xseclast_x>=0;
        xseclast_x   = xsec;

        // once the cross section stops increasing, I reduce the step size and
        // step backwards a little bit to handle cases that the max cross section
        // is grossly underestimated (very peaky distribution & large step)
        if(!increasing_x) {
          LOG("DISKinematics", pDEBUG) 
           << "d2xsec/dxdy|x stopped increasing. Stepping back & exiting x loop";
          double dlogxn = dlogx/(Nxb+1);
          //double dxn = dx/(Nxb+1);
          for(int ik=0; ik<Nxb; ik++) {
	     gx = TMath::Exp(TMath::Log(gx) - dlogxn);
   	     //gx = gx - dxn;
             interaction->KinePtr()->Setx(gx);
             kinematics::UpdateWQ2FromXY(interaction);

             xsec = fXSecModel->XSec(interaction, kPSxyfE);

             LOG("DISKinematics", pINFO) 
                << "xsec(y=" << gy << ", x=" << gx << ") = " << xsec;
	  }
          break;
        } // stepping back within last bin
     } // x
     increasing_y = max_xsec-xseclast_y>=0;
     xseclast_y   = max_xsec;
     if(!increasing_y) {
       LOG("DISKinematics", pDEBUG) 
                          << "d2xsec/dxdy stopped increasing. Exiting y loop";
       break;
     }
  }// y

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  //  max_xsec *= fSafetyFactor;
  //max_xsec *= ( (Ev<3.0) ? 2.5 : fSafetyFactor);
  max_xsec *= 3;

  SLOG("DISKinematics", pDEBUG) << interaction->AsString();
  SLOG("DISKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("DISKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________

