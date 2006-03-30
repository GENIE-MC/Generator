//____________________________________________________________________________
/*!

\class   genie::QELKinematicsGenerator

\brief   Generates values for the kinematic variables describing QEL neutrino
         interaction events.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Controls.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/QELKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::controls;

//___________________________________________________________________________
QELKinematicsGenerator::QELKinematicsGenerator() :
KineGeneratorWithCache("genie::QELKinematicsGenerator")
{

}
//___________________________________________________________________________
QELKinematicsGenerator::QELKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::QELKinematicsGenerator", config)
{

}
//___________________________________________________________________________
QELKinematicsGenerator::~QELKinematicsGenerator()
{

}
//___________________________________________________________________________
void QELKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Selects kinematic variables using the 'Rejection' method and adds them to
// the event record's summary

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Get the interaction and set the 'trust' bits
  Interaction * interaction = evrec->GetInteraction();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Get the limits for the generated Q2
  Range1D_t Q2 = this->Q2Range(interaction);

  if(Q2.max <=0 || Q2.min>=Q2.max) {
     LOG("QELKinematics", pWARN) << "No available phase space";
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

  //-- Try to select a valid Q2 using the rejection method
  register unsigned int iter = 0;

  assert(Q2.min>=0.);

  double e = 1E-6;
  double logQ2min = TMath::Log(Q2.min+e);
  double logQ2max = TMath::Log(Q2.max);
  double dlogQ2   = logQ2max - logQ2min;

  while(1) {

     iter++;
     if(iter > kRjMaxIterations) {
        LOG("QELKinematics", pWARN)
          << "Couldn't select a valid Q^2 after " << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kNoValidKinematics, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     // generate a Q2 value within the allowed phase space
     double gQ2 = TMath::Exp(logQ2min + dlogQ2 * rnd->Random1().Rndm());
     interaction->GetKinematicsPtr()->SetQ2(gQ2);

     LOG("QELKinematics", pINFO) << "Trying: Q^2 = " << gQ2;

     double xsec = fXSecModel->XSec(interaction);
     double t    = xsec_max * rnd->Random1().Rndm();

     LOG("QELKinematics", pINFO)
             << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert(xsec < xsec_max);

     if(t < xsec) {
        // kinematical selection done.
        LOG("QELKinematics", pINFO) << "Selected: Q^2 = " << gQ2;

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec);

        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);
        return;
     }
  }// iterations
}
//___________________________________________________________________________
void QELKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELKinematicsGenerator::LoadConfig(void)
{
// Load sub-algorithms and config data to reduce the number of registry
// lookups

  fXSecModel = dynamic_cast<const XSecAlgorithmI *> (
                            this->SubAlg("xsec-alg-name", "xsec-param-set"));
  assert(fXSecModel);

  //-- Get the user kinematical limits on Q2
  fQ2min = fConfig->GetDoubleDef("Q2-min", -999999);
  fQ2max = fConfig->GetDoubleDef("Q2-max",  999999);

  //-- Safety factor for the maximum differential cross section
  fSafetyFactor = fConfig->GetDoubleDef("max-xsec-safety-factor", 1.25);

  //-- Minimum energy for which max xsec would be cached, forcing explicit
  //   calculation for lower eneries
  fEMin = fConfig->GetDoubleDef("min-energy-cached", -1.0);
}
//____________________________________________________________________________
Range1D_t QELKinematicsGenerator::Q2Range(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction
  Range1D_t Q2 = utils::kinematics::Q2Range_M(interaction);
  LOG("QELKinematics", pDEBUG)
               << "Physical Q2 range = (" << Q2.min << ", " << Q2.max << ")";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  utils::kinematics::ApplyCutsToKineLimits(Q2, fQ2min, fQ2max);
  LOG("QELKinematics", pDEBUG)
      << "(Physical & User) Q2 range = (" << Q2.min << ", " << Q2.max << ")";

  return Q2;
}
//___________________________________________________________________________
double QELKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small dQ2 step.

  double max_xsec = 0.0;
  const int N = 40;

  const InitialState & init_state = interaction -> GetInitialState();
  double E = init_state.GetProbeE(kRfStruckNucAtRest);

  Range1D_t rQ2 = this->Q2Range(interaction);
  if( rQ2.max < 1e-3 || rQ2.min <=0 ) return 0.;
  if(E<0.6) utils::kinematics::ApplyCutsToKineLimits(rQ2, E/200., 1.2*E);

  const double logQ2min = TMath::Log(rQ2.min);
  const double logQ2max = TMath::Log(rQ2.max);
  const double dlogQ2   = (logQ2max - logQ2min) /(N-1);

  double xseclast = -1;
  bool   increasing;

  for(int i=0; i<N; i++) {
     double Q2 = TMath::Exp(logQ2min + i * dlogQ2);
     interaction->GetKinematicsPtr()->SetQ2(Q2);

     double xsec = fXSecModel->XSec(interaction);
     max_xsec = TMath::Max(xsec, max_xsec);

     increasing = xsec-xseclast>0;
     xseclast   = xsec;
     if(!increasing) break;
  }//Q^2

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  max_xsec *= fSafetyFactor;

  SLOG("QELKinematics", pDEBUG) << interaction->AsString();
  SLOG("QELKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("QELKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________

