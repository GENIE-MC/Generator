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

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

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

  Interaction * interaction = evrec->GetInteraction();

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  double xsec_max = this->MaxXSec(evrec);

  //-- Get the kinematical limits for the generated x,y
  Range1D_t y = this->yRange(interaction);
  const double logxmin = TMath::Log(kASmallNum);
  const double logxmax = TMath::Log(1.-kASmallNum);
  const double logymin = TMath::Log(y.min+kASmallNum);
  const double logymax = TMath::Log(y.max-kASmallNum);
  const double dlogx   = (logxmax - logxmin);
  const double dlogy   = (logymax - logymin);

  //------ Try to select a valid x,y pair
  register unsigned int iter = 0;

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
     double gx = TMath::Exp(logxmin + dlogx * rnd->Random1().Rndm());
     double gy = TMath::Exp(logymin + dlogy * rnd->Random1().Rndm());
     interaction->GetKinematicsPtr()->Setx(gx);
     interaction->GetKinematicsPtr()->Sety(gy);
     LOG("COHKinematics", pINFO)
                   << "Trying: (x = " << gx << ", y = " << gy << ")";

     double xsec    = fXSecModel->XSec(interaction, kPSlogxlogyfE);
     double tryxsec = xsec_max * rnd->Random1().Rndm();

     LOG("COHKinematics", pINFO)
          << "xsec: (computed) = " << xsec << ", (generated) = " << tryxsec;
     assert(xsec < xsec_max);

     if(tryxsec < xsec) {
        // kinematical selection done.
        LOG("COHKinematics", pINFO)
                             << "Selected: x = " << gx << ", y = " << gy;

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

        LOG("COHKinematics", pINFO)
          << "Selected: t = "<< gt << ", from ["<< tmin << ", "<< tmax << "]";

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

  const int N = 75;
  Range1D_t y = this->yRange(in);

  const double logxmin = TMath::Log(kASmallNum);
  const double logxmax = TMath::Log(1.-kASmallNum);
  const double logymin = TMath::Log(y.min+kASmallNum);
  const double logymax = TMath::Log(y.max-kASmallNum);
  const double dlogx   = (logxmax - logxmin) /(N-1);
  const double dlogy   = (logymax - logymin) /(N-1);

  double Ev  = in->GetInitialState().GetProbeE(kRfLab);

  for(int i=0; i<N; i++) {
   double gx = TMath::Exp(logxmin + i * dlogx);
   for(int j=0; j<N; j++) {
     double gy = TMath::Exp(logymin + j * dlogy);

     double Q2 = 2*kNucleonMass*gx*gy*Ev;
     if(Q2 > 0.5) continue;

     in->GetKinematicsPtr()->Setx(gx);
     in->GetKinematicsPtr()->Sety(gy);

     double xsec = fXSecModel->XSec(in, kPSlogxlogyfE);
     max_xsec = TMath::Max(max_xsec, xsec);
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

  fRo = fConfig->GetDoubleDef("Ro", gc->GetDouble("COH-Ro"));

  fSafetyFactor = fConfig->GetDoubleDef("max-xsec-safety-factor", 1.3);
  fEMin         = fConfig->GetDoubleDef("min-energy-cached",     -1.0);

  fXSecModel = dynamic_cast<const XSecAlgorithmI *> (
                            this->SubAlg("xsec-alg-name", "xsec-param-set"));
  assert(fXSecModel);
}
//____________________________________________________________________________

