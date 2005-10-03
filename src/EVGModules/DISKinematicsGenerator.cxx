//____________________________________________________________________________
/*!

\class   genie::DISKinematicsGenerator

\brief   Generates values for the kinematic variables describing DIS v
         interaction events.

         Is a concrete implementation of the EventRecordVisitorI interface.

         Part of its implementation, related with the caching and retrieval of
         previously computed values, is inherited from the KineGeneratorWithCache
         abstract class.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Constants.h"
#include "EVGModules/DISKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
DISKinematicsGenerator::DISKinematicsGenerator() :
KineGeneratorWithCache()
{
  fName = "genie::DISKinematicsGenerator";
}
//___________________________________________________________________________
DISKinematicsGenerator::DISKinematicsGenerator(const char * param_set) :
KineGeneratorWithCache(param_set)
{
  fName = "genie::DISKinematicsGenerator";

  this->FindConfig();
}
//___________________________________________________________________________
DISKinematicsGenerator::~DISKinematicsGenerator()
{

}
//___________________________________________________________________________
void DISKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Selects (x,y) kinematic variables using the 'Rejection' method and adds
// them to the event record's summary

  Interaction * interaction = evrec->GetInteraction();

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Get the selected cross section calculator
  const XSecAlgorithmI * xsec_alg =
            dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                                     "xsec-alg-name", "xsec-param-set"));

  //-- For the subsequent kinematic selection with the rejection method:
  //   Valculate the max differential cross section or retrieve it from
  //   the cache (if something similar was computed at a previous step).
  double xsec_max = this->MaxXSec(interaction);
  xsec_max *= 1.3;

  //------ Try to select a valid x,y pair

  //-- Get the limits for the generated (x,y) pairs
  Range1D_t x = this->XRange();
  Range1D_t y = this->YRange();

  register unsigned int iter = 0;

  while(1) {
     //-- generate ax x,y pair & set it to interaction
     double gx = x.min + (x.max-x.min) * rnd->Random2().Rndm();
     double gy = y.min + (y.max-y.min) * rnd->Random2().Rndm();

     interaction->GetScatParamsPtr()->Set("x", gx);
     interaction->GetScatParamsPtr()->Set("y", gy);

     LOG("DISKinematics", pINFO) << "Trying: x = "<< gx << ", y = "<< gy;

     //-- only proceed if I've got allowed kinematics
     if( this->ValidKinematics(interaction) ) {

         double xsec = xsec_alg->XSec(interaction);
         double t    = xsec_max * rnd->Random2().Rndm();

         LOG("DISKinematics", pINFO)
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
             LOG("DISKinematics", pERROR)
                << "*** Could not select a valid (x,y) pair after "
                                               << iter << " iterations";
             assert(false);
         }

     } // valid kinematics
  } // iterations
}
//___________________________________________________________________________
Range1D_t DISKinematicsGenerator::WRange(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction

  Range1D_t W = utils::kinematics::WRange(interaction);
  LOG("DISKinematics", pDEBUG)
       << "\n Physical W integration range: "
                                 << "[" << W.min << ", " << W.max << "] GeV";

  //-- Get the user kinematical limits
  double min = (fConfig->Exists("W-min")) ? fConfig->GetDouble("W-min") : -1;
  double max = (fConfig->Exists("W-max")) ? fConfig->GetDouble("W-max") : -1;

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if ( utils::math::IsWithinLimits(min, W) ) W.min = min;
  if ( utils::math::IsWithinLimits(max, W) ) W.max = max;

  LOG("DISKinematics", pDEBUG)
       << "\n (Physical && User) W integration range: "
                                 << "[" << W.min << ", " << W.max << "] GeV";
//  assert( W.min < W.max && W.min >= 0 );

  return W;
}
//___________________________________________________________________________
Range1D_t DISKinematicsGenerator::Q2Range(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction

  Range1D_t Q2 = utils::kinematics::Q2Range_xy(interaction);
  LOG("DISKinematics", pDEBUG)
       << "\n Physical Q2 integration range: "
                            << "[" << Q2.min << ", " << Q2.max << "] GeV^2";

  //-- Get the user kinematical limits
  double min = (fConfig->Exists("Q2-min")) ? fConfig->GetDouble("Q2-min") : -1;
  double max = (fConfig->Exists("Q2-max")) ? fConfig->GetDouble("Q2-max") : -1;

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if ( utils::math::IsWithinLimits(min, Q2) ) Q2.min = min;
  if ( utils::math::IsWithinLimits(max, Q2) ) Q2.max = max;

  LOG("DISKinematics", pDEBUG)
       << "\n (Physical && User) Q2 integration range: "
                            << "[" << Q2.min << ", " << Q2.max << "] GeV^2";
//  assert( Q2.min < Q2.max && Q2.min >= 0 );

  return Q2;
}
//___________________________________________________________________________
Range1D_t DISKinematicsGenerator::XRange(void) const
{
// Get the user-defined range for x or set the default range [0,1]

  Range1D_t x;

  x.min = fConfig->Exists("x-min") ? fConfig->GetDouble("x-min") : 0.;
  x.max = fConfig->Exists("x-max") ? fConfig->GetDouble("x-max") : 1.;

  assert(x.min >= 0 && x.max <= 1 && x.min < x.max);

  return x;
}
//___________________________________________________________________________
Range1D_t DISKinematicsGenerator::YRange(void) const
{
// Get the user-defined range for y or set the default range [0,1]

  Range1D_t y;

  y.min = fConfig->Exists("y-min") ? fConfig->GetDouble("y-min") : 0.;
  y.max = fConfig->Exists("y-max") ? fConfig->GetDouble("y-max") : 1.;

  assert(y.min >= 0 && y.max <= 1 && y.min < y.max);

  return y;
}
//___________________________________________________________________________
bool DISKinematicsGenerator::ValidKinematics(
                                       const Interaction * interaction) const
{
// Checks whether the selected kinematic variables x, y yield valid W, Q^2

  //-- get initial state information
  const InitialState & init_state = interaction->GetInitialState();
  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double Ev  = p4->Energy();
  double M   = init_state.GetTarget().StruckNucleonMass();
  double M2  = M*M;

  delete p4;

  //-- get current set of generated x,y
  double gx  = interaction->GetScatteringParams().x();
  double gy  = interaction->GetScatteringParams().y();

  //-- compute the corresponding generated W, Q^2
  double gW  = TMath::Sqrt( utils::math::NonNegative(M2+2*Ev*M*gy*(1-gx)) );
  double gQ2 = 2*gx*gy*M*Ev;

  //-- get the W, Q^2 range
  //   (the physical one, unless narrowed down by external user cuts)

  Range1D_t W  = this->WRange(interaction);
  Range1D_t Q2 = this->Q2Range(interaction);

  LOG("DISKinematics", pDEBUG)
                     << "Generated values: W = " << gW << ", Q^2 = " << gQ2;
  LOG("DISKinematics", pDEBUG)
                          << "W range = (" << W.min << ", " << W.max << ")";
  LOG("DISKinematics", pDEBUG)
                      << "Q^2 range = (" << Q2.min << ", " << Q2.max << ")";

  bool valid_W  = utils::math::IsWithinLimits( gW, W  );
  bool valid_Q2 = utils::math::IsWithinLimits( gQ2,Q2 );

  if(valid_W && valid_Q2) return true;
  else return false;
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

  const int Nx = 40;
  const int Ny = 40;

  double max_xsec = -1.0;

  Range1D_t rx = this->XRange();
  Range1D_t ry = this->YRange();

  const double logxmin = TMath::Log(rx.min);
  const double logxmax = TMath::Log(rx.max);
  const double logymin = TMath::Log(ry.min);
  const double logymax = TMath::Log(ry.max);
  const double dlogx   = (logxmax - logxmin) /(Nx-1);
  const double dlogy   = (logymax - logymin) /(Ny-1);

  const XSecAlgorithmI * xsec_alg =
            dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                                    "xsec-alg-name", "xsec-param-set"));
  for(int ix=0; ix<Nx; ix++) {

     double x = TMath::Exp(logxmin + ix * dlogx);
     interaction->GetScatParamsPtr()->Set("x", x);

     for(int iy=0; iy<Ny; iy++) {

         double y = TMath::Exp(logymin + iy * dlogy);
         interaction->GetScatParamsPtr()->Set("y", y);

         if( this->ValidKinematics(interaction) ) {

            double xsec = xsec_alg->XSec(interaction);
            max_xsec = TMath::Max(xsec, max_xsec);
         }
     } // y
  }// x

  SLOG("DISKinematics", pDEBUG) << interaction->AsString();
  SLOG("DISKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("DISKinematics", pDEBUG) << "Computed using alg = " << *xsec_alg;

  return max_xsec;
}
//___________________________________________________________________________

