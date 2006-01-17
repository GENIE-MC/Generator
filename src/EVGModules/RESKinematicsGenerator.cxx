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

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  double xsec_max = this->MaxXSec(evrec);

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
     //   (the physically allowed Q2's, unless an external cut is imposed)
     Range1D_t Q2 = this->Q2Range(interaction);
     if(Q2.min<=0. || Q2.min>Q2.max) continue;
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

     if(t < xsec) {
        // kinematical selection done.
        LOG("RESKinematics", pINFO)
                            << "Selected: W = " << gW << ", Q2 = " << gQ2;
        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec);
        // reset 'trust' bits
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
  Range1D_t Q2 = utils::kinematics::Q2Range_W(interaction);
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

  const int    NQ2 = 15;
  const double e   = 1e-4;
  const double MD  = 1.23;

  // Set W around the value where d^2xsec/dWdQ^2 peaks
  Range1D_t rW = this->WRange(interaction);
  const double W = (utils::math::IsWithinLimits(MD, rW)) ? MD : rW.max-e;
  interaction->GetKinematicsPtr()->SetW(W);

  // Set a Q2 range, within the allowed region (inclusing user cuts), in
  // which d^2xsec/dWdQ^2 peaks
  Range1D_t rQ2 = this->Q2Range(interaction);
  if( rQ2.max < kMinQ2Limit || rQ2.min <=0 ) return 0.;
  if(E<0.6) utils::kinematics::ApplyCutsToKineLimits(rQ2, 0.05*E, 1.5*E);

  const double logQ2min = TMath::Log(rQ2.min);
  const double logQ2max = TMath::Log(rQ2.max);
  const double dlogQ2   = (logQ2max - logQ2min) /(NQ2-1);
  for(int iq2=0; iq2<NQ2; iq2++) {
     double Q2 = TMath::Exp(logQ2min + iq2 * dlogQ2);
     interaction->GetKinematicsPtr()->SetQ2(Q2);
     double xsec = fXSecModel->XSec(interaction);
     max_xsec = TMath::Max(xsec, max_xsec);
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

