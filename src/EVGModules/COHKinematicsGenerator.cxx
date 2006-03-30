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
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/COHKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
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
void COHKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Selects kinematic variables using the 'Rejection' method and adds them to
// the event record's summary

  Interaction * interaction = evrec->GetInteraction();

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  double xsec_max = this->MaxXSec(evrec);

  //------ Try to select a valid x,y pair
  const double e = 1E-6;
  register unsigned int iter = 0;

  //-- Get the kinematical limits for the generated x,y
  Range1D_t x;
  x.min=0;
  x.max=1;
  Range1D_t y = this->yRange(interaction);

  const double logxmin = TMath::Log(x.min+e);
  const double logxmax = TMath::Log(x.max-e);
  const double rlogx   = (logxmax - logxmin);
  const double logymin = TMath::Log(y.min+e);
  const double logymax = TMath::Log(y.max-e);
  const double rlogy   = (logymax - logymin);

  while(1) {

     iter++;
     if(iter > kRjMaxIterations) {
        LOG("COHKinematics", pWARN)
             << "*** Could not select a valid (x,y) pair after "
                                               << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kNoValidKinematics, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     double gx = TMath::Exp(logxmin + rlogx * rnd->Random1().Rndm());
     double gy = TMath::Exp(logymin + rlogy * rnd->Random1().Rndm());
     interaction->GetKinematicsPtr()->Setx(gx);
     interaction->GetKinematicsPtr()->Sety(gy);
     LOG("COHKinematics", pINFO)
                   << "Trying: (x = " << gx << ", y = " << gy << ")";

     double xsec = fXSecModel->XSec(interaction);
     double t    = xsec_max * rnd->Random1().Rndm();

     LOG("COHKinematics", pINFO)
             << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert(xsec < xsec_max);

     if(t < xsec) {
        // kinematical selection done.
        LOG("COHKinematics", pINFO)
                             << "Selected: x = " << gx << ", y = " << gy;
        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec);
        return;
     }
  }// iterations
}
//___________________________________________________________________________
Range1D_t COHKinematicsGenerator::yRange(const Interaction * in) const
{
  double Ev  = in->GetInitialState().GetProbeE(kRfLab);
  double Mpi = kPionMass;

  Range1D_t y;
  y.min = Mpi/Ev;
  y.max = 1.;

  LOG("COHKinematics", pDEBUG)
                  << "Physical y range = (" << y.min << ", " << y.max << ")";
 return y;
}
//___________________________________________________________________________
double COHKinematicsGenerator::ComputeMaxXSec(const Interaction * in) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.

  SLOG("COHKinematics", pDEBUG)
          << "Scanning the allowed phase space {K} for the max(dxsec/d{K})";

  double max_xsec = 0.;

  const double e = 1E-3;
  const int    N = 50;

  Range1D_t x;
  x.min=0.;
  x.max=1.;
  Range1D_t y = this->yRange(in);

  const double logxmin = TMath::Log(x.min+e);
  const double logxmax = TMath::Log(x.max-e);
  const double dlogx   = (logxmax - logxmin) /(N-1);
  const double logymin = TMath::Log(y.min+e);
  const double logymax = TMath::Log(y.max-e);
  const double dlogy   = (logymax - logymin) /(N-1);

  double Ev  = in->GetInitialState().GetProbeE(kRfLab);

  for(int i=0; i<N; i++) {
   double gx = TMath::Exp(logxmin + i * dlogx);
   for(int j=0; j<N; j++) {
     double gy = TMath::Exp(logymin + j * dlogy);

     double Q2 = 2*kNucleonMass*gx*gy*Ev;
     if(Q2 > 0.5) continue;

     in->GetKinematicsPtr()->Setx(gx);
     in->GetKinematicsPtr()->Sety(gy);

     double xsec = fXSecModel->XSec(in);
     max_xsec = TMath::Max(max_xsec, xsec);
   }//y
  }//x

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

  SLOG("COHKinematics", pDEBUG) << in->AsString();
  SLOG("COHKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("COHKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________
double COHKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->GetInitialState();
  double E = init_state.GetProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void COHKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void COHKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void COHKinematicsGenerator::LoadSubAlg(void)
{
  fXSecModel = dynamic_cast<const XSecAlgorithmI *> (
                            this->SubAlg("xsec-alg-name", "xsec-param-set"));
  assert(fXSecModel);
}
//____________________________________________________________________________
void COHKinematicsGenerator::LoadConfigData(void)
{
  fSafetyFactor = fConfig->GetDoubleDef("max-xsec-safety-factor", 1.3);
}
//____________________________________________________________________________

