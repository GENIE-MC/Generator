//____________________________________________________________________________
/*!

\class   genie::RESKinematicsGenerator

\brief   Generates resonance event (v+N->l+Resonance) kinematics.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 18, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "BaryonResonance/BaryonResonance.h"
#include "Conventions/Controls.h"
#include "EVGModules/RESKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
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
void RESKinematicsGenerator::ProcessEventRecord(GHepRecord * event_rec) const
{
// Selects (Q^2,W) kinematic variables using the 'Rejection' method and adds
// them to the event record's summary

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Valculate the max differential cross section or retrieve it from the
  //   cache (if something similar was computed at a previous step).
  Interaction * interaction = event_rec->GetInteraction();

  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  double xsec_max = this->MaxXSec(interaction);
  xsec_max *= fSafetyFactor;

  //-- Try to select a valid W, Q2 pair
  register unsigned int iter = 0;
  double e = 1E-6;

  //-- Compute the W limits
  //  (the physically allowed W's, unless an external cut is imposed)
  Range1D_t W = this->WRange(interaction);
  assert(W.min>0.);
  double logWmin  = TMath::Log(W.min+e);
  double logWmax  = TMath::Log(W.max);
  double dlogW    = logWmax - logWmin;

  while(1) {
     //-- Get a random W within its allowed limits
     double gW = TMath::Exp(logWmin + dlogW  * rnd->Random2().Rndm());
     interaction->GetKinematicsPtr()->SetW(gW);

     //-- Compute the allowed Q^2 limits for the selected W
     //   (the physically allowed W's, unless an external cut is imposed)
     Range1D_t Q2 = this->Q2Range(interaction);
     assert(Q2.min>0.);
     double logQ2min = TMath::Log(Q2.min+e);
     double logQ2max = TMath::Log(Q2.max);
     double dlogQ2   = logQ2max - logQ2min;

     //-- Get a random Q2 within its allowed limits
     //double gQ2 = Q2.min + (Q2.max - Q2.min) * rnd->Random2().Rndm();
     double gQ2 = TMath::Exp(logQ2min + dlogQ2 * rnd->Random2().Rndm());
     interaction->GetKinematicsPtr()->SetQ2(gQ2);

     LOG("RESKinematics", pINFO) << "Trying: W = " << gW << ", Q2 = " << gQ2;

     //--  The cross section algorithm that I am about to run should
     //    compute d^2xsec /dW dQ^2 for the *list* of currently considered
     //    baryon resonances and returns their sum weighted with the value
     //    of their Breit-Wigner distribution at the current W.
     double xsec = fXSecModel->XSec(interaction);
     double t    = xsec_max * rnd->Random2().Rndm();

     LOG("RESKinematics", pINFO)
              << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert(xsec < xsec_max);

     if( t < xsec ) {
        // kinematical selection done.
        LOG("RESKinematics", pINFO)
                            << "Selected: W = " << gW << ", Q2 = " << gQ2;
        // set the cross section for the selected kinematics
        event_rec->SetDiffXSec(xsec);

        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);
        return;
     }

     iter++;
     if(iter > kRjMaxIterations) {
         LOG("RESKinematics", pFATAL)
              << "*** Could not select a valid (W,Q^2) pair after "
                                                    << iter << " iterations";
         abort();
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
  fWmin = fConfig->GetDoubleDef("W-min", -1);
  fWmax = fConfig->GetDoubleDef("W-max", -1);

  //-- Get the user kinematical limits on Q2
  fQ2min = fConfig->GetDoubleDef("Q2-min", -1);
  fQ2max = fConfig->GetDoubleDef("Q2-max", -1);

  //-- Safety factor for the maximum differential cross section
  fSafetyFactor = fConfig->GetDoubleDef("max-xsec-safety-factor", 1.25);
}
//____________________________________________________________________________
Range1D_t RESKinematicsGenerator::WRange(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction
  Range1D_t W = utils::kinematics::WRange(interaction);
  LOG("RESKinematics", pDEBUG)
       << "\n Physical W integration range: "
                                 << "[" << W.min << ", " << W.max << "] GeV";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if ( utils::math::IsWithinLimits(fWmin, W) ) W.min = fWmin;
  if ( utils::math::IsWithinLimits(fWmax, W) ) W.max = fWmax;

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
  Range1D_t Q2 = utils::kinematics::Q2Range_W(interaction);
  LOG("RESKinematics", pDEBUG)
       << "\n Physical Q2 integration range: "
                            << "[" << Q2.min << ", " << Q2.max << "] GeV^2";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if ( utils::math::IsWithinLimits(fQ2min, Q2) ) Q2.min = fQ2min;
  if ( utils::math::IsWithinLimits(fQ2max, Q2) ) Q2.max = fQ2max;

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

  //const int NW  = 10;
  const int NQ2 = 20;

  double max_xsec = -1.0;

  //Range1D_t rW = this->WRange(interaction);

  //const double dW  = (rW.max-rW.min)/(NW-1);

  //for(int iw=0; iw<NW; iw++) {
  //   double W = rW.min + iw * dW;
     double W=1.232;
     interaction->GetKinematicsPtr()->SetW(W);

     Range1D_t rQ2 = this->Q2Range(interaction);
     const double logQ2min = TMath::Log(rQ2.min);
     const double logQ2max = TMath::Log(rQ2.max);
     const double dlogQ2   = (logQ2max - logQ2min) /(NQ2-1);

     for(int iq2=0; iq2<NQ2; iq2++) {
        double Q2 = TMath::Exp(logQ2min + iq2 * dlogQ2);
        interaction->GetKinematicsPtr()->SetQ2(Q2);

        double xsec = fXSecModel->XSec(interaction);

        max_xsec = TMath::Max(xsec, max_xsec);
     } // Q2
  //}// W

  LOG("RESKinematics", pDEBUG) << interaction->AsString();
  LOG("RESKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  LOG("RESKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________

