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

#include "Conventions/Controls.h"
#include "EVGModules/QELKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
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
void QELKinematicsGenerator::ProcessEventRecord(GHepRecord * event_rec) const
{
// Selects (x,y) kinematic variables using the 'Rejection' method and adds
// them to the event record's summary

  Interaction * interaction = event_rec->GetInteraction();

  //-- Get the selected cross section calculator
  const XSecAlgorithmI * xsec_alg =
               dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                                        "xsec-alg-name", "xsec-param-set"));
  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Valculate the max differential cross section or retrieve it from
  //   the cache (if something similar was computed at a previous step).
  double xsec_max = this->MaxXSec(interaction);
  xsec_max *= 1.3;

  //------ Try to select a valid Q2

  //-- Get the limits for the generated Q2

  Range1D_t Q2 = this->Q2Range(interaction);

  register unsigned int iter = 0;
  while(1) {
     double gQ2 = Q2.min + (Q2.max-Q2.min) * rnd->Random2().Rndm();
     interaction->GetKinematicsPtr()->SetQ2(gQ2);

     LOG("QELKinematics", pINFO) << "Trying: Q^2 = " << gQ2;

     double xsec = xsec_alg->XSec(interaction);
     double t    = xsec_max * rnd->Random2().Rndm();

     LOG("QELKinematics", pINFO)
             << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert( xsec < xsec_max );

     if( t < xsec ) {
        // kinematical selection done.
        LOG("QELKinematics", pINFO) << "Selected: Q^2 = " << gQ2;

        // set the cross section for the selected kinematics
        event_rec->SetDiffXSec(xsec);
        return;
     }

     iter++;
     if(iter > kRjMaxIterations) {
        LOG("QELKinematics", pERROR)
                  << "*** Could not select a valid (x,y) pair after "
                                                    << iter << " iterations";
        assert(false);
     }
  }// iterations
}
//___________________________________________________________________________
Range1D_t QELKinematicsGenerator::Q2Range(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction

  Range1D_t Q2 = utils::kinematics::Q2Range_M(interaction);

  LOG("QELKinematics", pDEBUG)
               << "Physical Q2 range = (" << Q2.min << ", " << Q2.max << ")";

  //-- Get the user kinematical limits
  double min = (fConfig->Exists("Q2-min")) ? fConfig->GetDouble("Q2-min") : -1;
  double max = (fConfig->Exists("Q2-max")) ? fConfig->GetDouble("Q2-max") : -1;

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if ( utils::math::IsWithinLimits(min, Q2) ) Q2.min = min;
  if ( utils::math::IsWithinLimits(max, Q2) ) Q2.max = max;

  LOG("QELKinematics", pDEBUG)
      << "(Physical & User) Q2 range = (" << Q2.min << ", " << Q2.max << ")";

//  assert( Q2.min < Q2.max && Q2.min >= 0 );

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

  const int N = 20;
  double max_xsec = -1.0;

  const XSecAlgorithmI * xsec_alg =
           dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                                     "xsec-alg-name", "xsec-param-set"));

  Range1D_t rQ2 = this->Q2Range(interaction);
  const double logQ2min = TMath::Log(rQ2.min);
  const double logQ2max = TMath::Log(rQ2.max);
  const double dlogQ2   = (logQ2max - logQ2min) /(N-1);

  for(int i=0; i<N; i++) {

     double Q2 = TMath::Exp(logQ2min + i * dlogQ2);
     interaction->GetKinematicsPtr()->SetQ2(Q2);

     double xsec = xsec_alg->XSec(interaction);

     max_xsec = TMath::Max(xsec, max_xsec);
  }//Q^2

  SLOG("QELKinematics", pDEBUG) << interaction->AsString();
  SLOG("QELKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("QELKinematics", pDEBUG) << "Computed using alg = " << *xsec_alg;

  return max_xsec;
}
//___________________________________________________________________________

