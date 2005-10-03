//____________________________________________________________________________
/*!

\class   genie::IMDKinematicsGenerator

\brief   Generates values for the kinematic variables describing inverse muon
         decay (IMD) neutrino interaction events.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 13, 2005

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "EVGModules/IMDKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/MathUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
IMDKinematicsGenerator::IMDKinematicsGenerator() :
KineGeneratorWithCache()
{
  fName = "genie::IMDKinematicsGenerator";
}
//___________________________________________________________________________
IMDKinematicsGenerator::IMDKinematicsGenerator(const char * param_set) :
KineGeneratorWithCache(param_set)
{
  fName = "genie::IMDKinematicsGenerator";

  this->FindConfig();
}
//___________________________________________________________________________
IMDKinematicsGenerator::~IMDKinematicsGenerator()
{

}
//___________________________________________________________________________
void IMDKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Selects (x,y) kinematic variables using the 'Rejection' method and adds
// them to the event record's summary

  Interaction * interaction = evrec->GetInteraction();

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

  //------ Try to select a valid inelastisity y

  double ymin = 0.0; // the xsec algorithm would internally compute the
  double ymax = 1.0; // kinematically allowd range, and return 0 if outside
  double dy   = ymax-ymin;

  register unsigned int iter = 0;
  while(1) {
     double y = ymin + dy * rnd->Random2().Rndm();
     interaction->GetScatParamsPtr()->Set("y", y);

     LOG("IMDKinematics", pINFO) << "Trying: y = " << y;

     double xsec = xsec_alg->XSec(interaction);
     double t    = xsec_max * rnd->Random2().Rndm();

     LOG("IMDKinematics", pINFO)
           << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert( xsec < xsec_max );

     if( t < xsec ) {
        // kinematical selection done.
        LOG("IMDKinematics", pINFO) << "Selected: y = " << y;

        // set the cross section for the selected kinematics
        interaction->SetDiffXSec(xsec);
        return;
     }

     iter++;
     if(iter > kRjMaxIterations) {
        LOG("IMDKinematics", pERROR)
              << "*** Could not select a valid y after "
                                            << iter << " iterations";
        assert(false);
     }
  }// iterations
}
//___________________________________________________________________________
double IMDKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small y step.

  const int N = 20;
  double max_xsec = -1.0;

  const XSecAlgorithmI * xsec_alg =
           dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                                     "xsec-alg-name", "xsec-param-set"));

  const double ymin = 0.;
  const double ymax = 1.;
  const double dy   = (ymax-ymin)/(N-1);

  for(int i=0; i<N; i++) {

     double y = ymin + i * dy;
     interaction->GetScatParamsPtr()->Set("y", y);

     double xsec = xsec_alg->XSec(interaction);

     max_xsec = TMath::Max(xsec, max_xsec);
  }//y

  SLOG("IMDKinematics", pDEBUG) << interaction->AsString();
  SLOG("IMDKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("IMDKinematics", pDEBUG) << "Computed using alg = " << *xsec_alg;

  return max_xsec;
}
//___________________________________________________________________________

