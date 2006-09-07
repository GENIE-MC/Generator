//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Controls.h"
#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EVGThreadException.h"
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

  //-- Get the interaction 
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);

  //-- Get neutrino energy and hit 'nucleon mass' 
  const InitialState & init_state = interaction->InitState();
  double Ev  = init_state.ProbeE(kRfHitNucRest);
  double M   = init_state.Tgt().HitNucP4().M(); // can be off m-shell

  //-- Get the physical W range taking into account any user cuts
  Range1D_t W  = this->WRange(interaction);
  if(W.max <=0 || W.min>=W.max) {
     LOG("DISKinematics", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kNoAvailablePhaseSpace, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("No available phase space");
     exception.SwitchOnFastForward();
     throw exception;
  }

  //-- Get the minimum Q^2 (corresponding to W = Wmin)
  interaction->KinePtr()->SetW(W.min);
  Range1D_t Q2  = this->Q2Range(interaction);

  //-- Set the corresponding x,y limits
  fXmax = 1. - kASmallNum;
  fYmax = 1. - kASmallNum;
  fXmin = 0;
  fYmin = 0;
  kinematics::WQ2toXY(Ev,M,W.min,Q2.min,fXmin, fYmin);

  LOG("DISKinematics", pINFO) << "x: [" << fXmin << ", " << fXmax << "]";
  LOG("DISKinematics", pINFO) << "y: [" << fYmin << ", " << fYmax << "]";

  assert(fXmin>0 && fYmin>0);

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //-- Try to select a valid (x,y) pair using the rejection method

  double dx = fXmax - fXmin;
  double dy = fYmax - fYmin;
  /*
  double logxmin = TMath::Log(fXmin);
  double logxmax = TMath::Log(fXmax);
  double logymin = TMath::Log(fYmin);
  double logymax = TMath::Log(fYmax);
  double dlogx   = logxmax-logxmin;
  double dlogy   = logymax-logymin;
  */
  double gx=-1, gy=-1, gW=-1, gQ2=-1, xsec=-1;

  register unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
       LOG("DISKinematics", pWARN)
        << " Couldn't select kinematics after " << iter << " iterations";

       evrec->EventFlags()->SetBitNumber(kNoValidKinematics, true);
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select kinematics");
       exception.SwitchOnFastForward();
       throw exception;
     }

     //-- random x,y
     gx = fXmin  + dx  * rnd->RndKine().Rndm();
     gy = fYmin  + dy  * rnd->RndKine().Rndm();

     LOG("DISKinematics", pINFO) 
                           << "Trying: x = " << gx << ", y = " << gy;

     interaction->KinePtr()->Setx(gx);
     interaction->KinePtr()->Sety(gy);

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
         	             << "Selected:  x = " << gx << ", y = " << gy;

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
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void DISKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void DISKinematicsGenerator::LoadSubAlg(void)
{
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed
  fXSecModel = dynamic_cast<const XSecAlgorithmI *> (
                            this->SubAlg("xsec-alg-name", "xsec-param-set"));
  assert(fXSecModel);
}
//____________________________________________________________________________
void DISKinematicsGenerator::LoadConfigData(void)
{
// Reads its configuration data from its configuration Registry and loads them
// in private data members to avoid looking up at the Registry all the time.

  //-- Get the user kinematical limits on W
  fWminCut = fConfig->GetDoubleDef("W-min", -999999);
  fWmaxCut = fConfig->GetDoubleDef("W-max",  999999);

  //-- Get the user kinematical limits on Q2
  fQ2minCut = fConfig->GetDoubleDef("Q2-min", -999999);
  fQ2maxCut = fConfig->GetDoubleDef("Q2-max",  999999);

  //-- Safety factor for the maximum differential cross section
  fSafetyFactor = fConfig->GetDoubleDef("max-xsec-safety-factor", 1.25);

  //-- Minimum energy for which max xsec would be cached, forcing explicit
  //   calculation for lower eneries
  fEMin = fConfig->GetDoubleDef("min-energy-cached", 1.0);

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  fMaxXSecDiffTolerance = fConfig->GetDoubleDef("max-xsec-diff-tolerance",0.);
  assert(fMaxXSecDiffTolerance>=0);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  fGenerateUniformly = fConfig->GetBoolDef("uniform-over-phase-space", false);
}
//____________________________________________________________________________
Range1D_t DISKinematicsGenerator::WRange(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction
  Range1D_t W = kinematics::KineRange(interaction,kKVW);
  LOG("DISKinematics", pDEBUG)
             << "Physical W range: " << "[" << W.min << ", " << W.max << "]";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if(fWminCut>0 && fWmaxCut>0) 
         kinematics::ApplyCutsToKineLimits(W, fWminCut, fWmaxCut);
  LOG("DISKinematics", pDEBUG)
    << "W range (including cuts): " << "[" << W.min << ", " << W.max << "]";

  return W;
}
//___________________________________________________________________________
Range1D_t DISKinematicsGenerator::Q2Range(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction
  Range1D_t Q2 = kinematics::KineRange(interaction, kKVQ2);
  LOG("DISKinematics", pDEBUG)
         << "Physical Q2 range: " << "[" << Q2.min << ", " << Q2.max << "]";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if(fQ2minCut>0 && fQ2maxCut>0)
         kinematics::ApplyCutsToKineLimits(Q2, fQ2minCut, fQ2maxCut);
  LOG("DISKinematics", pDEBUG)
     << "Q2 range (including cuts): " << "["<< Q2.min<< ", "<< Q2.max<< "]";

  return Q2;
}
//___________________________________________________________________________
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
    ywindow = .5;
    Ny  = 20;
    Nx  = 80;
    Nxb = 10;
  }

  if(tgt.HitQrkIsSet() && tgt.HitSeaQrk()) {
    xpeak    = .1;
    xwindow  = .1;
  }

  double xmin    = TMath::Max(xpeak-xwindow, fXmin);
  double xmax    = TMath::Min(xpeak+xwindow, fXmax);
  double ymin    = TMath::Max(ypeak-ywindow, fYmin);
  double ymax    = TMath::Min(ypeak+ywindow, fYmax);
  double logxmin = TMath::Log10(xmin);
  double logxmax = TMath::Log10(xmax);
  double logymin = TMath::Log10(ymin);
  double logymax = TMath::Log10(ymax);
  double dlogx   = (logxmax-logxmin)/(Nx-1);
  double dlogy   = (logymax-logymin)/(Ny-1);

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
  max_xsec *= ( (Ev<3.5) ? 2. : fSafetyFactor);

  SLOG("DISKinematics", pDEBUG) << interaction->AsString();
  SLOG("DISKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("DISKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________

