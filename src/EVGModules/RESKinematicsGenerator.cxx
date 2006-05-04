//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - November 18, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Controls.h"
#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/RESKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"

using namespace genie;
using namespace genie::controls;

//___________________________________________________________________________
RESKinematicsGenerator::RESKinematicsGenerator() :
KineGeneratorWithCache("genie::RESKinematicsGenerator")
{

}
//___________________________________________________________________________
RESKinematicsGenerator::RESKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::RESKinematicsGenerator", config)
{

}
//___________________________________________________________________________
RESKinematicsGenerator::~RESKinematicsGenerator()
{

}
//___________________________________________________________________________
void RESKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Selects kinematic variables using the 'Rejection' method and adds them to
// the event record's summary

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Get the interaction and set the 'trust' bits
  Interaction * interaction = evrec->GetInteraction();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Compute the W limits
  //  (the physically allowed W's, unless an external cut is imposed)
  Range1D_t W = this->WRange(interaction);
  assert(W.min>=0. && W.min<W.max);

  if(W.max <=0 || W.min>=W.max) {
     LOG("RESKinematics", pWARN) << "No available phase space";
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

  //-- Try to select a valid W, Q2 pair using the rejection method
  register unsigned int iter = 0;

  double Wmin = W.min + kASmallNum;
  double Wmax = W.max - kASmallNum;
  double dW   = Wmax - Wmin;

  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
         LOG("RESKinematics", pWARN)
              << "*** Could not select a valid (W,Q^2) pair after "
                                                    << iter << " iterations";
         evrec->EventFlags()->SetBitNumber(kNoValidKinematics, true);
         genie::exceptions::EVGThreadException exception;
         exception.SetReason("Couldn't select kinematics");
         exception.SwitchOnFastForward();
         throw exception;
     }

     //-- Get a random W within its allowed limits
     double gW = Wmin + dW  * rnd->Random1().Rndm();
     interaction->GetKinematicsPtr()->SetW(gW);

     //-- Compute the allowed Q^2 limits for the selected W
     //   (the physically allowed Q2's, unless an external cut is imposed)
     Range1D_t Q2 = this->Q2Range(interaction);
     if(Q2.max<=0. || Q2.min>=Q2.max) continue;
     double logQ2min = TMath::Log(Q2.min+kASmallNum);
     double logQ2max = TMath::Log(Q2.max-kASmallNum);
     double dlogQ2   = logQ2max - logQ2min;

     //-- Get a random Q2 within its allowed limits
     double gQ2 = TMath::Exp(logQ2min + dlogQ2 * rnd->Random1().Rndm());
     interaction->GetKinematicsPtr()->SetQ2(gQ2);

     LOG("RESKinematics", pINFO) << "Trying: W = " << gW << ", Q2 = " << gQ2;

     //--  The cross section algorithm that I am about to run should
     //    compute d^2xsec /dW dQ^2 for the *list* of currently considered
     //    baryon resonances and returns their sum weighted with the value
     //    of their Breit-Wigner distribution at the current W.
     double xsec = fXSecModel->XSec(interaction, kPSWlogQ2fE);
     double t    = xsec_max * rnd->Random1().Rndm();

     LOG("RESKinematics", pINFO)
              << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert(xsec < xsec_max);

     if(t < xsec) {
        // kinematical selection done.
        LOG("RESKinematics", pINFO)
                            << "Selected: W = " << gW << ", Q2 = " << gQ2;

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec);

        // reset 'trust' bits
        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);

        // compute x,y for selected W,Q2
        // note: hit nucleon can be off the mass-shell
        double gx=-1, gy=-1;
        const InitialState & init_state = interaction->GetInitialState();
        double E = init_state.GetProbeE(kRfStruckNucAtRest);
        double M = init_state.GetTarget().StruckNucleonP4()->M(); 
        utils::kinematics::WQ2toXY(E,M,gW,gQ2,gx,gy);

        // lock selected kinematics & clear running values
        interaction->GetKinematicsPtr()->SetQ2(gQ2, true);
        interaction->GetKinematicsPtr()->SetW (gW,  true);
        interaction->GetKinematicsPtr()->Setx (gx,  true);
        interaction->GetKinematicsPtr()->Sety (gy,  true);
        interaction->GetKinematicsPtr()->ClearRunningValues();

        return;
     }
  } // iterations
}
//___________________________________________________________________________
void RESKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void RESKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void RESKinematicsGenerator::LoadSubAlg(void)
{
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed
  fXSecModel = dynamic_cast<const XSecAlgorithmI *> (
                            this->SubAlg("xsec-alg-name", "xsec-param-set"));
  assert(fXSecModel);
}
//____________________________________________________________________________
void RESKinematicsGenerator::LoadConfigData(void)
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
Range1D_t RESKinematicsGenerator::WRange(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction
  Range1D_t W = utils::kinematics::KineRange(interaction, kKVW);
  LOG("RESKinematics", pDEBUG)
       << "\n Physical W integration range: "
                                 << "[" << W.min << ", " << W.max << "] GeV";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  utils::kinematics::ApplyCutsToKineLimits(W, fWmin, fWmax);
  LOG("RESKinematics", pDEBUG)
       << "\n (Physical & User) W integration range: "
                                 << "[" << W.min << ", " << W.max << "] GeV";
  return W;
}
//___________________________________________________________________________
Range1D_t RESKinematicsGenerator::Q2Range(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction
  Range1D_t Q2 = utils::kinematics::KineRange(interaction, kKVQ2);
  LOG("RESKinematics", pDEBUG)
       << "\n Physical Q2 integration range: "
                            << "[" << Q2.min << ", " << Q2.max << "] GeV^2";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  utils::kinematics::ApplyCutsToKineLimits(Q2, fQ2min, fQ2max);
  LOG("RESKinematics", pDEBUG)
       << "\n (Physical && User) Q2 integration range: "
                            << "[" << Q2.min << ", " << Q2.max << "] GeV^2";
  return Q2;
}
//___________________________________________________________________________
double RESKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But this needs to be fast - do not use a very fine grid.

  double max_xsec = 0.;

  const InitialState & init_state = interaction -> GetInitialState();
  double E = init_state.GetProbeE(kRfStruckNucAtRest);

  double md;
  if(!interaction->GetExclusiveTag().KnownResonance()) md=1.23;
  else {
    Resonance_t res = interaction->GetExclusiveTag().Resonance();
    md=utils::res::Mass(res);
  }

  const int    NQ2 = 15;
  const double MD  = md;

  // Set W around the value where d^2xsec/dWdQ^2 peaks
  Range1D_t rW = this->WRange(interaction);
  double W;
  if(utils::math::IsWithinLimits(MD, rW))  W = MD;
  else {
    if (MD>=rW.max) W = rW.max-kASmallNum;
    else            W = rW.min+kASmallNum;
  }
  interaction->GetKinematicsPtr()->SetW(W);

  // Set a Q2 range, within the allowed region (inclusing user cuts), in
  // which d^2xsec/dWdQ^2 peaks
  Range1D_t rQ2 = this->Q2Range(interaction);
  if( rQ2.max < kMinQ2Limit || rQ2.min <=0 ) return 0.;

  const double logQ2min = TMath::Log(rQ2.min);
  const double logQ2max = TMath::Log(rQ2.max);
  const double dlogQ2   = (logQ2max - logQ2min) /(NQ2-1);

  double xseclast = -1;
  bool   increasing;

  for(int iq2=0; iq2<NQ2; iq2++) {
     double Q2 = TMath::Exp(logQ2min + iq2 * dlogQ2);
     interaction->GetKinematicsPtr()->SetQ2(Q2);
     double xsec = fXSecModel->XSec(interaction, kPSWlogQ2fE);
     max_xsec = TMath::Max(xsec, max_xsec);

     increasing = xsec-xseclast>0;
     xseclast=xsec;
     if(!increasing) break;
  } // Q2

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  // Apply larger safety factor for smaller energies.
  max_xsec *= ( (E<0.6) ? 2. : fSafetyFactor);

  LOG("RESKinematics", pDEBUG) << interaction->AsString();
  LOG("RESKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  LOG("RESKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________

