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
KineGeneratorWithCache()
{
  fName = "genie::RESKinematicsGenerator";
}
//___________________________________________________________________________
RESKinematicsGenerator::RESKinematicsGenerator(const char * param_set) :
KineGeneratorWithCache(param_set)
{
  fName = "genie::RESKinematicsGenerator";

  this->FindConfig();
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

  //-- Get the selected cross section calculator
  const XSecAlgorithmI * xsec_alg =
             dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                                         "xsec-alg-name", "xsec-param-set"));

  //-- For the subsequent kinematic selection with the rejection method:
  //   Valculate the max differential cross section or retrieve it from the
  //   cache (if something similar was computed at a previous step).

  Interaction * interaction = event_rec->GetInteraction();

  double xsec_max = this->MaxXSec(interaction);
  xsec_max *= 1.3;

  //-- Try to select a valid W, Q2 pair

  //-- Compute the W limits
  //  (the physically allowed W's, unless an external cut is imposed)

  Range1D_t W = this->WRange(interaction);

  register unsigned int iter = 0;

  while(1) {
     //-- Get a random W within its allowed limits

     double gW = W.min + (W.max - W.min) * rnd->Random2().Rndm();
     interaction->GetScatParamsPtr()->Set("W", gW);

     //-- Compute the allowed Q^2 limits for the selected W
     //   (the physically allowed W's, unless an external cut is imposed)

     Range1D_t Q2 = this->Q2Range(interaction);

     //-- Get a random Q2 within its allowed limits

     double gQ2 = Q2.min + (Q2.max - Q2.min) * rnd->Random2().Rndm();
     interaction->GetScatParamsPtr()->Set("Q2", gQ2);

     LOG("RESKinematics", pINFO) << "Trying: W = " << gW << ", Q2 = " << gQ2;

     //--  The cross section algorithm that I am about to run should
     //    compute d^2xsec /dW dQ^2 for the *list* of currently considered
     //    baryon resonances and returns their sum weighted with the value
     //    of their Breit-Wigner distribution at the current W.

     double xsec = xsec_alg->XSec(interaction);
     double t    = xsec_max * rnd->Random2().Rndm();

     LOG("RESKinematics", pINFO)
              << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert(xsec < xsec_max);

     if( t < xsec ) {
        // kinematical selection done.
        LOG("RESKinematics", pINFO)
                            << "Selected: W = " << gW << ", Q2 = " << gQ2;
        // set the cross section for the selected kinematics
        interaction->SetDiffXSec(xsec);
        return;
     }

     iter++;
     if(iter > kRjMaxIterations) {
         LOG("RESKinematics", pERROR)
              << "*** Could not select a valid (W,Q^2) pair after "
                                                    << iter << " iterations";
         assert(false);
     }
  } // iterations
}
//___________________________________________________________________________
Range1D_t RESKinematicsGenerator::WRange(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction

  Range1D_t W = utils::kinematics::WRange(interaction);
  LOG("RESKinematics", pDEBUG)
       << "\n Physical W integration range: "
                                 << "[" << W.min << ", " << W.max << "] GeV";

  //-- Get the user kinematical limits
  double min = (fConfig->Exists("W-min")) ? fConfig->GetDouble("W-min") : -1;
  double max = (fConfig->Exists("W-max")) ? fConfig->GetDouble("W-max") : -1;

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.

  if ( utils::math::IsWithinLimits(min, W) ) W.min = min;
  if ( utils::math::IsWithinLimits(max, W) ) W.max = max;

  LOG("RESKinematics", pDEBUG)
       << "\n (Physical & User) W integration range: "
                                 << "[" << W.min << ", " << W.max << "] GeV";
//  assert( W.min < W.max && W.min >= 0 );

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

  //-- Get the user kinematical limits
  double min = (fConfig->Exists("Q2-min")) ? fConfig->GetDouble("Q2-min") : -1;
  double max = (fConfig->Exists("Q2-max")) ? fConfig->GetDouble("Q2-max") : -1;

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.

  if ( utils::math::IsWithinLimits(min, Q2) ) Q2.min = min;
  if ( utils::math::IsWithinLimits(max, Q2) ) Q2.max = max;

  LOG("RESKinematics", pDEBUG)
       << "\n (Physical && User) Q2 integration range: "
                            << "[" << Q2.min << ", " << Q2.max << "] GeV^2";
//  assert( Q2.min < Q2.max && Q2.min >= 0 );

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

  const int NW  = 10;
  const int NQ2 = 20;

  double max_xsec = -1.0;

  const XSecAlgorithmI * xsec_alg =
           dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                                       "xsec-alg-name", "xsec-param-set"));

  Range1D_t rW = this->WRange(interaction);

  const double dW  = (rW.max-rW.min)/(NW-1);

  for(int iw=0; iw<NW; iw++) {
     double W = rW.min + iw * dW;
     interaction->GetScatParamsPtr()->Set("W", W);

     Range1D_t rQ2 = this->Q2Range(interaction);
     const double logQ2min = TMath::Log(rQ2.min);
     const double logQ2max = TMath::Log(rQ2.max);
     const double dlogQ2   = (logQ2max - logQ2min) /(NQ2-1);

     for(int iq2=0; iq2<NQ2; iq2++) {
        double Q2 = TMath::Exp(logQ2min + iq2 * dlogQ2);
        interaction->GetScatParamsPtr()->Set("Q2", Q2);

        double xsec = xsec_alg->XSec(interaction);

        max_xsec = TMath::Max(xsec, max_xsec);
     } // Q2
  }// W

  LOG("RESKinematics", pDEBUG) << interaction->AsString();
  LOG("RESKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  LOG("RESKinematics", pDEBUG) << "Computed using alg = " << *xsec_alg;

  return max_xsec;
}
//___________________________________________________________________________

