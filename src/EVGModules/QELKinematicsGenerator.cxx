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

  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Valculate the max differential cross section or retrieve it from
  //   the cache (if something similar was computed at a previous step).
  double xsec_max = this->MaxXSec(interaction);
  xsec_max *= fSafetyFactor;

  //------ Try to select a valid Q2

  //-- Get the limits for the generated Q2
  Range1D_t Q2 = this->Q2Range(interaction);

  register unsigned int iter = 0;
  while(1) {
     double gQ2 = Q2.min + (Q2.max-Q2.min) * rnd->Random2().Rndm();
     interaction->GetKinematicsPtr()->SetQ2(gQ2);

     LOG("QELKinematics", pINFO) << "Trying: Q^2 = " << gQ2;

     double xsec = fXSecModel->XSec(interaction);
     double t    = xsec_max * rnd->Random2().Rndm();

     LOG("QELKinematics", pINFO)
             << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert( xsec < xsec_max );

     if( t < xsec ) {
        // kinematical selection done.
        LOG("QELKinematics", pINFO) << "Selected: Q^2 = " << gQ2;

        // set the cross section for the selected kinematics
        event_rec->SetDiffXSec(xsec);

        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);
        return;
     }

     iter++;
     if(iter > kRjMaxIterations) {
        LOG("QELKinematics", pFATAL)
                  << "*** Could not select a valid (x,y) pair after "
                                                    << iter << " iterations";
        abort();
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
  fQ2min = fConfig->GetDoubleDef("Q2-min", -1);
  fQ2max = fConfig->GetDoubleDef("Q2-max", -1);

  //-- Safety factor for the maximum differential cross section
  fSafetyFactor = fConfig->GetDoubleDef("max-xsec-safety-factor", 1.25);
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
  if ( utils::math::IsWithinLimits(fQ2min, Q2) ) Q2.min = fQ2min;
  if ( utils::math::IsWithinLimits(fQ2max, Q2) ) Q2.max = fQ2max;

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

  const int N = 40;
  double max_xsec = -1.0;

  Range1D_t rQ2 = this->Q2Range(interaction);
  const double logQ2min = TMath::Log(rQ2.min);
  const double logQ2max = TMath::Log(rQ2.max);
  const double dlogQ2   = (logQ2max - logQ2min) /(N-1);

  for(int i=0; i<N; i++) {

     double Q2 = TMath::Exp(logQ2min + i * dlogQ2);
     interaction->GetKinematicsPtr()->SetQ2(Q2);

     double xsec = fXSecModel->XSec(interaction);

     max_xsec = TMath::Max(xsec, max_xsec);
  }//Q^2

  SLOG("QELKinematics", pDEBUG) << interaction->AsString();
  SLOG("QELKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("QELKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________

