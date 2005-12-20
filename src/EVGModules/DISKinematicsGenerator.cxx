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
#include "Conventions/Controls.h"
#include "EVGModules/DISKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::controls;

//___________________________________________________________________________
DISKinematicsGenerator::DISKinematicsGenerator() :
KineGeneratorWithCache("genie::DISKinematicsGenerator")
{

}
//___________________________________________________________________________
DISKinematicsGenerator::DISKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::DISKinematicsGenerator", config)
{

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

  //-- For the subsequent kinematic selection with the rejection method:
  //   Valculate the max differential cross section or retrieve it from
  //   the cache (if something similar was computed at a previous step).
  double xsec_max = this->MaxXSec(interaction);
  xsec_max *= 1.3;

  //------ Try to select a valid x,y pair
  register unsigned int iter = 0;
  while(1) {
     //-- generate ax x,y pair & set it to interaction
     double gx = fXmin + fdX * rnd->Random2().Rndm();
     double gy = fYmin + fdY * rnd->Random2().Rndm();

     interaction->GetKinematicsPtr()->Setx(gx);
     interaction->GetKinematicsPtr()->Sety(gy);

     LOG("DISKinematics", pINFO) << "Trying: x = "<< gx << ", y = "<< gy;

     //-- only proceed if I've got allowed kinematics
     if( this->ValidKinematics(interaction) ) {

         double xsec = fXSecModel->XSec(interaction);
         double t    = xsec_max * rnd->Random2().Rndm();

         LOG("DISKinematics", pINFO)
             << "xsec: (computed) = " << xsec << ", (generated) = " << t;
         assert( xsec < xsec_max );

         if( t < xsec ) {
            // kinematical selection done.
            LOG("DISKinematics", pINFO)
                             << "Selected: x = " << gx << ", y = " << gy;
            // set the cross section for the selected kinematics
            evrec->SetDiffXSec(xsec);
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
void DISKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void DISKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void DISKinematicsGenerator::LoadSubAlg(void)
{
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed
  fXSecModel = dynamic_cast<const XSecAlgorithmI *> (
                            this->SubAlg("xsec-alg-name", "xsec-param-set"));
  assert(fXSecModel);
}
//____________________________________________________________________________
void DISKinematicsGenerator::LoadConfigData(void)
{
// Reads its configuration data from its configuration Registry and loads them
// in private data members to avoid looking up at the Registry all the time.

  //-- Get the user kinematical limits on W
  fWmin = fConfig->GetDoubleDef("W-min", -1);
  fWmax = fConfig->GetDoubleDef("W-max", -1);

  //-- Get the user kinematical limits on Q2
  fQ2min = fConfig->GetDoubleDef("Q2-min", -1);
  fQ2max = fConfig->GetDoubleDef("Q2-max", -1);

  double t0 = 1E-6;
  double t1 = 1. - t0;

  fXmin = fConfig->GetDoubleDef("x-min", t0);
  fXmax = fConfig->GetDoubleDef("x-max", t1);
  fYmin = fConfig->GetDoubleDef("y-min", t0);
  fYmax = fConfig->GetDoubleDef("y-max", t1);

  assert(fXmin > 0 && fXmax < 1 && fXmin < fXmax);
  assert(fYmin > 0 && fYmax < 1 && fYmin < fYmax);

  fdX = fXmax - fXmin;
  fdY = fYmax - fYmin;
}
//____________________________________________________________________________
Range1D_t DISKinematicsGenerator::WRange(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction
  Range1D_t W = utils::kinematics::WRange(interaction);
  LOG("DISKinematics", pDEBUG)
       << "\n Physical W integration range: "
                                 << "[" << W.min << ", " << W.max << "] GeV";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if ( utils::math::IsWithinLimits(fWmin, W) ) W.min = fWmin;
  if ( utils::math::IsWithinLimits(fWmax, W) ) W.max = fWmax;

  LOG("DISKinematics", pDEBUG)
       << "\n (Physical && User) W integration range: "
                                 << "[" << W.min << ", " << W.max << "] GeV";
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

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if ( utils::math::IsWithinLimits(fQ2min, Q2) ) Q2.min = fQ2min;
  if ( utils::math::IsWithinLimits(fQ2max, Q2) ) Q2.max = fQ2max;

  LOG("DISKinematics", pDEBUG)
       << "\n (Physical && User) Q2 integration range: "
                            << "[" << Q2.min << ", " << Q2.max << "] GeV^2";
  return Q2;
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
  double gx  = interaction->GetKinematics().x();
  double gy  = interaction->GetKinematics().y();

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

  const double logxmin = TMath::Log(fXmin);
  const double logxmax = TMath::Log(fXmax);
  const double logymin = TMath::Log(fYmin);
  const double logymax = TMath::Log(fYmax);
  const double dlogx   = (logxmax - logxmin) /(Nx-1);
  const double dlogy   = (logymax - logymin) /(Ny-1);

  for(int ix=0; ix<Nx; ix++) {

     double x = TMath::Exp(logxmin + ix * dlogx);
     interaction->GetKinematicsPtr()->Setx(x);

     for(int iy=0; iy<Ny; iy++) {
         double y = TMath::Exp(logymin + iy * dlogy);
         interaction->GetKinematicsPtr()->Sety(y);

         if( this->ValidKinematics(interaction) ) {

            double xsec = fXSecModel->XSec(interaction);
            max_xsec = TMath::Max(xsec, max_xsec);
         }
     } // y
  }// x

  SLOG("DISKinematics", pDEBUG) << interaction->AsString();
  SLOG("DISKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("DISKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________

