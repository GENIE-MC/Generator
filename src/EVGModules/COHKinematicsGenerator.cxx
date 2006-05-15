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

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/COHKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
COHKinematicsGenerator::COHKinematicsGenerator() :
KineGeneratorWithCache("genie::COHKinematicsGenerator")
{

}
//___________________________________________________________________________
COHKinematicsGenerator::COHKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::COHKinematicsGenerator", config)
{

}
//___________________________________________________________________________
COHKinematicsGenerator::~COHKinematicsGenerator()
{

}
//___________________________________________________________________________
void COHKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Selects kinematic variables using the 'Rejection' method and adds them to
// the event record's summary

  if(fGenerateUniformly) {
    LOG("COHKinematics", pNOTICE)
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

  //-- Get the kinematical limits for the generated x,y
  Interaction * interaction = evrec->GetInteraction();
  Range1D_t y = this->yRange(interaction);

  const double xmin = kASmallNum;
  const double xmax = 1.- kASmallNum;
  const double ymin = y.min + kASmallNum;
  const double ymax = y.max - kASmallNum;
  const double dx   = xmax - xmin;
  const double dy   = ymax - ymin;

  //------ Try to select a valid x,y pair

  register unsigned int iter = 0;
  bool accept=false;
  double xsec=-1;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("COHKinematics", pWARN)
             << "*** Could not select a valid (x,y) pair after "
                                               << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kNoValidKinematics, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }
     double gx = xmin + dx * rnd->Random1().Rndm();
     double gy = ymin + dy * rnd->Random1().Rndm();
     interaction->GetKinematicsPtr()->Setx(gx);
     interaction->GetKinematicsPtr()->Sety(gy);
     LOG("COHKinematics", pINFO) << "Trying: x = " << gx << ", y = " << gy;

     // computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSxyfE);

     //-- decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);

        double t = xsec_max * rnd->Random1().Rndm();
        LOG("COHKinematics", pINFO)
                              << "xsec= " << xsec << ", J= 1, Rnd= " << t;
        accept = (t<xsec);
     }
     else { 
        accept = (xsec>0);
     }

     if(accept) {
         // ----------------- KINEMATICAL SELECTION DONE -------------------

        LOG("COHKinematics", pNOTICE) << "Selected: x = "<< gx << ", y = "<< gy;

        // the COH cross section should be a triple differential cross section
        // d^2xsec/dxdydt where t is the the square of the 4p transfer to the
        // nucleus. The cross section used for kinematical selection should have
        // the t-dependence integrated out. The t-dependence is of the form
        // ~exp(-bt). Now that the x,y kinematical variables have been selected
        // we can generate a t using the t-dependence as a PDF.
        const InitialState & init_state = interaction->GetInitialState();
        double Ev    = init_state.GetProbeE(kRfLab);
        double Epi   = gy*Ev; // pion energy
        double Epi2  = TMath::Power(Epi,2);
        double pme2  = kPionMass2/Epi2;   
        double xME   = kNucleonMass*gx/Epi;
        double tA    = 1. + xME - 0.5*pme2;
        double tB    = TMath::Sqrt(1.+ 2*xME) * TMath::Sqrt(1.-pme2);
        double tmin  = 2*Epi2 * (tA-tB);
        double tmax  = 2*Epi2 * (tA+tB);
        double A     = (double) init_state.GetTarget().A(); 
        double A13   = TMath::Power(A,1./3.);
        double R     = fRo * A13; // nuclear radius
        double R2    = TMath::Power(R,2.);
        double b     = 0.33333 * R2;
        double tsum  = (TMath::Exp(-b*tmin) - TMath::Exp(-b*tmax))/b; 
        double rt    = tsum * rnd->Random1().Rndm();
        double gt    = -1.*TMath::Log(-1.*b*rt + TMath::Exp(-1.*b*tmin))/b;

        LOG("COHKinematics", pNOTICE)
          << "Selected: t = "<< gt << ", from ["<< tmin << ", "<< tmax << "]";

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
        if(fGenerateUniformly) {
          double vol     = y.max-y.min; // dx=1, dt: irrelevant
          double totxsec = evrec->GetXSec();
          double wght    = (vol/totxsec)*xsec;
          LOG("COHKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

          // apply computed weight to the current event weight
          wght *= evrec->GetWeight();
          LOG("COHKinematics", pNOTICE) << "Current event wght = " << wght;
          evrec->SetWeight(wght);
        }

        // lock selected kinematics & clear running values
        interaction->GetKinematicsPtr()->Setx(gx, true);
        interaction->GetKinematicsPtr()->Sety(gy, true);
        interaction->GetKinematicsPtr()->Sett(gt, true);
        interaction->GetKinematicsPtr()->SetW(kPionMass, true);
        interaction->GetKinematicsPtr()->SetQ2(2*kNucleonMass*gx*gy*Ev, true);
        interaction->GetKinematicsPtr()->ClearRunningValues();

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec);

        return;
     }
  }// iterations
}
//___________________________________________________________________________
Range1D_t COHKinematicsGenerator::yRange(const Interaction * in) const
{
  double Ev  = in->GetInitialState().GetProbeE(kRfLab);
  double Mpi = kPionMass;

  Range1D_t y;
  y.min = Mpi/Ev;
  y.max = 1.;

  LOG("COHKinematics", pDEBUG)
                  << "Physical y range = (" << y.min << ", " << y.max << ")";
 return y;
}
//___________________________________________________________________________
double COHKinematicsGenerator::ComputeMaxXSec(const Interaction * in) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.

  SLOG("COHKinematics", pDEBUG)
          << "Scanning the allowed phase space {K} for the max(dxsec/d{K})";

  double max_xsec = 0.;

  const int Nx  =  80;
  const int Ny  =  25;
  const int Nyb =  15;

  Range1D_t y = this->yRange(in);

  const double logxmin = TMath::Log(kASmallNum);
  const double logxmax = TMath::Log(1.-kASmallNum);
  const double logymin = TMath::Log(y.min+kASmallNum);
  const double logymax = TMath::Log(y.max-kASmallNum);
  const double dlogx   = (logxmax - logxmin) /(Nx-1);
  const double dlogy   = (logymax - logymin) /(Ny-1);

  double Ev  = in->GetInitialState().GetProbeE(kRfLab);

  double xseclast = -1;
  bool   increasing;

  for(int i=0; i<Nx; i++) {
   double gx = TMath::Exp(logxmin + i * dlogx);
   for(int j=0; j<Ny; j++) {
     double gy = TMath::Exp(logymin + j * dlogy);

     double Q2 = 2*kNucleonMass*gx*gy*Ev;
     if(Q2 > 0.1) continue;

     in->GetKinematicsPtr()->Setx(gx);
     in->GetKinematicsPtr()->Sety(gy);

     double xsec = fXSecModel->XSec(in, kPSxyfE);
     LOG("COHKinematics", pDEBUG)  
     	              << "xsec(x= " << gx << ", y= " << gy << ") = " << xsec;

     max_xsec = TMath::Max(max_xsec, xsec);
     increasing = xsec-xseclast>=0;
     xseclast   = xsec;

     // once the xsec stops increasing, I reduce the step size and
     // step backwards a little bit to handle cases that the max xsec
     // is grossly underestimated (very peaky distribution & large step)
     if(!increasing) {
       double dlogyb = dlogy/(Nyb+1);
       for(int jb=0; jb<Nyb; jb++) {
	 gy = TMath::Exp(TMath::Log(gy) - dlogyb);
         if(gy < y.min) continue;
         in->GetKinematicsPtr()->Sety(gy);
         xsec = fXSecModel->XSec(in, kPSxyfE);
         LOG("COHKinematics", pDEBUG)  
     	              << "xsec(x= " << gx << ", y= " << gy << ") = " << xsec;
         max_xsec = TMath::Max(xsec, max_xsec);
       }
       break;
     }
   }//y
  }//x

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

  SLOG("COHKinematics", pDEBUG) << in->AsString();
  SLOG("COHKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("COHKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________
double COHKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->GetInitialState();
  double E = init_state.GetProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void COHKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHKinematicsGenerator::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  //-- COH model parameter Ro
  fRo = fConfig->GetDoubleDef("Ro", gc->GetDouble("COH-Ro"));

  //-- max xsec safety factor (for rejection method) and min cached energy
  fSafetyFactor = fConfig->GetDoubleDef("max-xsec-safety-factor", 1.2);
  fEMin         = fConfig->GetDoubleDef("min-energy-cached",     -1.0);

  //-- Differential cross section model
  fXSecModel = dynamic_cast<const XSecAlgorithmI *> (
                            this->SubAlg("xsec-alg-name", "xsec-param-set"));

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  fGenerateUniformly = fConfig->GetBoolDef("uniform-over-phase-space", false);

  assert(fXSecModel);
}
//____________________________________________________________________________

