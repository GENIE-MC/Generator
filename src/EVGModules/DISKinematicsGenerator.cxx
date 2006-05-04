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

using namespace genie;
using namespace genie::controls;

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
// Selects kinematic variables using the 'Rejection' method and adds them to
// the event record's summary

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Get the interaction and set the 'trust' bits
  Interaction * interaction = evrec->GetInteraction();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

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

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  double xsec_max = this->MaxXSec(evrec);

  //-- get energy and 'nucleon mass' needed for W,Q2->x,y conversion
  const InitialState & init_state = interaction->GetInitialState();
  double Ev  = init_state.GetProbeE(kRfStruckNucAtRest);
  double M   = init_state.GetTarget().StruckNucleonP4()->M(); // can be off m-shell

  //-- Try to select a valid (x,y) pair using the rejection method
  register unsigned int iter = 0;

  double logxmin  = TMath::Log(kMinX);
  double logxmax  = TMath::Log(kMaxX);
  double logymin  = TMath::Log(kMinY);
  double logymax  = TMath::Log(kMaxY);
  double dlogx    = logxmax - logxmin;
  double dlogy    = logymax - logymin;

  double gx=-1, gy=-1, gW=-1, gQ2=-1;

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
     gx = TMath::Exp(logxmin  + dlogx  * rnd->Random1().Rndm());
     gy = TMath::Exp(logymin  + dlogy  * rnd->Random1().Rndm());

     //-- (x,y) => (W,Q2)
     utils::kinematics::XYtoWQ2(Ev,M,gW,gQ2,gx,gy);

     interaction->GetKinematicsPtr()->Setx(gx);
     interaction->GetKinematicsPtr()->Sety(gy);
     interaction->GetKinematicsPtr()->SetW(gW);
     interaction->GetKinematicsPtr()->SetQ2(gQ2);

     LOG("DISKinematics", pINFO) 
             << "Trying: (x=" << gx << ", y=" << gy 
                          << ") -> (W=" << gW << ", Q2=" << gQ2 << ")";

     //-- check whether we are in the allowed phase space
     Range1D_t Q2 = this->Q2Range(interaction);
     bool inW  = utils::math::IsWithinLimits(gW, W );
     bool inQ2 = utils::math::IsWithinLimits(gQ2,Q2);
     bool in   = inW && inQ2;

     if(!in) {
       LOG("DISKinematics", pINFO) << "*** Not an allowed phase space point";
       continue;
     }

     //-- compute the cross section for current kinematics
     double xsec = fXSecModel->XSec(interaction, kPSxyfE);

     //-- accept current kinematics?
     double t = xsec_max * rnd->Random1().Rndm();
     LOG("DISKinematics", pINFO)
             << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert(xsec < xsec_max);
     if(t < xsec) {
         // kinematical selection done.
         LOG("DISKinematics", pNOTICE)
                  << "Selected: W = "<< gW << ", Q2 = "<< gQ2
                           << " (=> x = " << gx << ", y = " << gy << ")";

         interaction->ResetBit(kISkipProcessChk);
         interaction->ResetBit(kISkipKinematicChk);

         // set the cross section for the selected kinematics
         evrec->SetDiffXSec(xsec);

         // lock selected kinematics & clear running values
         interaction->GetKinematicsPtr()->SetW (gW,  true);
         interaction->GetKinematicsPtr()->SetQ2(gQ2, true);
         interaction->GetKinematicsPtr()->Setx (gx,  true);
         interaction->GetKinematicsPtr()->Sety (gy,  true);
         interaction->GetKinematicsPtr()->ClearRunningValues();

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
  fWmin = fConfig->GetDoubleDef("W-min", -999999);
  fWmax = fConfig->GetDoubleDef("W-max",  999999);

  //-- Get the user kinematical limits on Q2
  fQ2min = fConfig->GetDoubleDef("Q2-min", -999999);
  fQ2max = fConfig->GetDoubleDef("Q2-max",  999999);

  //-- Safety factor for the maximum differential cross section
  fSafetyFactor = fConfig->GetDoubleDef("max-xsec-safety-factor", 1.25);

  //-- Minimum energy for which max xsec would be cached, forcing explicit
  //   calculation for lower eneries
  fEMin = fConfig->GetDoubleDef("min-energy-cached", 1.0);

}
//____________________________________________________________________________
Range1D_t DISKinematicsGenerator::WRange(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction
  Range1D_t W = utils::kinematics::KineRange(interaction,kKVW);
  LOG("DISKinematics", pDEBUG)
       << "\n Physical W integration range: "
                                 << "[" << W.min << ", " << W.max << "] GeV";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if(fWmin>0 && fWmax>0)
      utils::kinematics::ApplyCutsToKineLimits(W, fWmin, fWmax);
  LOG("DISKinematics", pDEBUG)
       << "\n (Physical && User) W integration range: "
                                 << "[" << W.min << ", " << W.max << "] GeV";
  return W;
}
//___________________________________________________________________________
Range1D_t DISKinematicsGenerator::Q2Range(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction
  Range1D_t Q2 = utils::kinematics::KineRange(interaction, kKVQ2);
  LOG("DISKinematics", pDEBUG)
       << "\n Physical Q2 integration range: "
                            << "[" << Q2.min << ", " << Q2.max << "] GeV^2";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if(fQ2min>0 && fQ2max>0)
      utils::kinematics::ApplyCutsToKineLimits(Q2, fQ2min, fQ2max);
  LOG("DISKinematics", pDEBUG)
       << "\n (Physical && User) Q2 integration range: "
                            << "[" << Q2.min << ", " << Q2.max << "] GeV^2";
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

  double max_xsec = 0.0;

  const InitialState & init_state = interaction -> GetInitialState();
  double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
  double M  = init_state.GetTarget().StruckNucleonMass(); // on m-shell here

  const int    Nx       = 30;
  const int    Ny       = 20;
  const double logxmin  = TMath::Log(kMinX);
  const double logxmax  = TMath::Log(kMaxX);
  const double logymin  = TMath::Log(kMinY);
  const double logymax  = TMath::Log(kMaxY);
  const double dlogx    = (logxmax - logxmin)/Nx;
  const double dlogy    = (logymax - logymin)/Ny;

  LOG("DISKinematics", pDEBUG)
                       << "Computing max xsec in allowed x,y phase space";

  //-- Get the physical W range taking into account any user cuts
  Range1D_t W = this->WRange(interaction);
  LOG("DISKinematics", pNOTICE)
                       << "W range = (" << W.min << ", " << W.max << ")";
  assert(W.min>0 && W.max>W.min);

  double gW=-1, gQ2=-1, gx=-1, gy=-1;

  for(int ix=0; ix<Nx; ix++) {
     gx = TMath::Exp(logxmin + ix*dlogx);

     double xseclast = -1;
     bool increasing;

     for(int iy=0; iy<Ny; iy++) {
        gy = TMath::Exp(logymin + iy*dlogy);

        //-- (x,y) => (W,Q2)
        utils::kinematics::XYtoWQ2(Ev,M,gW,gQ2,gx,gy);

        interaction->GetKinematicsPtr()->Setx(gx);
        interaction->GetKinematicsPtr()->Sety(gy);
        interaction->GetKinematicsPtr()->SetW(gW);
        interaction->GetKinematicsPtr()->SetQ2(gQ2);

        //-- check whether we are in the allowed phase space
        Range1D_t Q2 = this->Q2Range(interaction);
        bool inW  = utils::math::IsWithinLimits(gW, W );
        bool inQ2 = utils::math::IsWithinLimits(gQ2,Q2);
        bool in   = inW && inQ2;

        if(!in) continue;

        // update maximum xsec
        double xsec = fXSecModel->XSec(interaction, kPSxyfE);
        max_xsec = TMath::Max(xsec, max_xsec);

        increasing = xsec-xseclast>0;
        xseclast   = xsec;
	//        if(!increasing) break;
     } // y
  }// x

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  max_xsec *= fSafetyFactor;

  SLOG("DISKinematics", pDEBUG) << interaction->AsString();
  SLOG("DISKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("DISKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________

