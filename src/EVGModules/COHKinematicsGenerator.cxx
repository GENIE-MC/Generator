//____________________________________________________________________________
/*!

\class   genie::COHKinematicsGenerator

\brief   Generates values for the kinematic variables describing QEL neutrino
         interaction events.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "EVGModules/COHKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
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
void COHKinematicsGenerator::ProcessEventRecord(GHepRecord * event_rec) const
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

  //------ Try to select a valid x,y pair

  //-- Get the kinematical limits for the generated x,y

  Range1D_t x;
  x.min=0.;
  x.max=1.;

  Range1D_t y = this->yRange(interaction);

  register unsigned int iter = 0;
  while(1) {
     double gx = x.min + (x.max-x.min) * rnd->Random2().Rndm();
     double gy = y.min + (y.max-y.min) * rnd->Random2().Rndm();
     interaction->GetScatParamsPtr()->Set("x", gx);
     interaction->GetScatParamsPtr()->Set("y", gy);
     LOG("COHKinematics", pINFO)
                       << "Trying: (x = " << gx << ", y = " << gy << ")";

     double xsec = xsec_alg->XSec(interaction);
     double t    = xsec_max * rnd->Random2().Rndm();

     LOG("COHKinematics", pINFO)
             << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert( xsec < xsec_max );

     if( t < xsec ) {
        // kinematical selection done.
        LOG("DISKinematics", pINFO)
                             << "Selected: x = " << gx << ", y = " << gy;

        // set the cross section for the selected kinematics
        interaction->SetDiffXSec(xsec);
        return;
     }

     iter++;
     if(iter > kRjMaxIterations) {
        LOG("COHKinematics", pERROR)
             << "*** Could not select a valid (x,y) pair after "
                                               << iter << " iterations";
        assert(false);
     }
  }// iterations
}
//___________________________________________________________________________
Range1D_t COHKinematicsGenerator::yRange(
                                      const Interaction * interaction) const
{
  Range1D_t y;

  double Ev  = interaction->GetInitialState().GetProbeE(kRfLab);
  double Mpi = kPionMass;

  y.min = Mpi/Ev;
  y.max = 1.;

  LOG("COHKinematics", pDEBUG)
                  << "Physical y range = (" << y.min << ", " << y.max << ")";

 return y;
}
//___________________________________________________________________________
double COHKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.

  const int N = 20;
  double max_xsec = -1.0;

  const XSecAlgorithmI * xsec_alg =
           dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                                      "xsec-alg-name", "xsec-param-set"));
  Range1D_t x;
  x.min=0.;
  x.max=1.;

  Range1D_t y = this->yRange(interaction);

  const double logxmin = TMath::Log(x.min);
  const double logxmax = TMath::Log(x.max);
  const double dlogx   = (logxmax - logxmin) /(N-1);

  const double logymin = TMath::Log(y.min);
  const double logymax = TMath::Log(y.max);
  const double dlogy   = (logymax - logymin) /(N-1);

  for(int i=0; i<N; i++) {
   double gx = TMath::Exp(logxmin + i * dlogx);
   interaction->GetScatParamsPtr()->Set("x", gx);
   for(int j=0; j<N; j++) {
     double gy = TMath::Exp(logymin + j * dlogy);
     interaction->GetScatParamsPtr()->Set("y", gy);

     max_xsec = TMath::Max(max_xsec, xsec_alg->XSec(interaction));
   }//y
  }//x

  SLOG("COHKinematics", pDEBUG) << interaction->AsString();
  SLOG("COHKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("COHKinematics", pDEBUG) << "Computed using alg = " << *xsec_alg;

  return max_xsec;
}
//___________________________________________________________________________

