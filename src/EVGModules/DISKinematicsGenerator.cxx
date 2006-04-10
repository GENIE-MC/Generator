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

#include <TMath.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Controls.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/DISKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
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
// Selects kinematic variables using the 'Rejection' method and adds them to
// the event record's summary

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Get the interaction and set the 'trust' bits
  Interaction * interaction = evrec->GetInteraction();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Get the physical W range taking into account any user cuts
  Range1D_t W  = this->WRange(interaction);

  if(W.max <=0 || W.min>=W.max) {
     LOG("DISKinematics", pWARN) << "No available phase space";
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

  //-- Try to select a valid W,Q2 (=>x,y) pair using the rejection method
  register unsigned int iter = 0;
  double e = 1E-6;

  assert(W.min>=0.);
  double logWmin  = TMath::Log(W.min+e);
  double logWmax  = TMath::Log(W.max);
  double dlogW    = logWmax - logWmin;

  while(1) {

     iter++;
     if(iter > kRjMaxIterations) {
       LOG("DISKinematics", pWARN)
        << " Couldn't select kinematics after " << iter << " iterations";

       evrec->EventFlags()->SetBitNumber(kNoValidKinematics, true);
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select kinematics");
       exception.SwitchOnFastForward();
       throw exception;
     }

     //-- generate a W value within the allowed phase space
     double gW  = TMath::Exp(logWmin  + dlogW  * rnd->Random1().Rndm());
     interaction->GetKinematicsPtr()->SetW(gW);

     //-- Get the physical Q2 range (for current W) taking into account
     //   any user cuts
     Range1D_t Q2 = this->Q2Range(interaction);
     if(Q2.max <=0 || Q2.min>Q2.max) continue;
     double logQ2min = TMath::Log(Q2.min+e);
     double logQ2max = TMath::Log(Q2.max);
     double dlogQ2   = logQ2max - logQ2min;

     //-- generate a Q2 value within the allowed phase space
     double gQ2 = TMath::Exp(logQ2min + dlogQ2 * rnd->Random1().Rndm());
     interaction->GetKinematicsPtr()->SetQ2(gQ2);

     LOG("DISKinematics", pINFO) << "Trying: W = "<< gW << ", Q2 = "<< gQ2;

     //-- W,Q2 => x,y
     this->SetKineXY(interaction);
     double gx = interaction->GetKinematics().x();
     double gy = interaction->GetKinematics().y();

     //-- compute the cross section for current kinematics
     double xsec = fXSecModel->XSec(interaction);

     //-- accept current kinematics?
     double t = xsec_max * rnd->Random1().Rndm();
     LOG("DISKinematics", pINFO)
             << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert(xsec < xsec_max);
     if(t < xsec) {
         // kinematical selection done.
         LOG("DISKinematics", pINFO)
                  << "Selected: W = "<< gW << ", Q2 = "<< gQ2
                           << " (=> x = " << gx << ", y = " << gy << ")";

         interaction->ResetBit(kISkipProcessChk);
         interaction->ResetBit(kISkipKinematicChk);

         // set the cross section for the selected kinematics
         evrec->SetDiffXSec(xsec);

         // lock selected kinematics & clear running values
         interaction->GetKinematicsPtr()->SetW (gW,  true);
         interaction->GetKinematicsPtr()->SetQ2(gQ2, true);
         interaction->GetKinematicsPtr()->Setx (gx,  true);
         interaction->GetKinematicsPtr()->Sety (gy,  true);
         interaction->GetKinematicsPtr()->ClearRunningValues();

         return;
     }
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
  fWmin = fConfig->GetDoubleDef("W-min", -999999);
  fWmax = fConfig->GetDoubleDef("W-max",  999999);

  //-- Get the user kinematical limits on Q2
  fQ2min = fConfig->GetDoubleDef("Q2-min", -999999);
  fQ2max = fConfig->GetDoubleDef("Q2-max",  999999);

  //-- Safety factor for the maximum differential cross section
  fSafetyFactor = fConfig->GetDoubleDef("max-xsec-safety-factor", 1.25);

  //-- Minimum energy for which max xsec would be cached, forcing explicit
  //   calculation for lower eneries
  fEMin = fConfig->GetDoubleDef("min-energy-cached", 1.0);

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
  if(fWmin>0 && fWmax>0)
      utils::kinematics::ApplyCutsToKineLimits(W, fWmin, fWmax);
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
  Range1D_t Q2 = utils::kinematics::Q2Range_W(interaction);
  LOG("DISKinematics", pDEBUG)
       << "\n Physical Q2 integration range: "
                            << "[" << Q2.min << ", " << Q2.max << "] GeV^2";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if(fQ2min>0 && fQ2max>0)
      utils::kinematics::ApplyCutsToKineLimits(Q2, fQ2min, fQ2max);
  LOG("DISKinematics", pDEBUG)
       << "\n (Physical && User) Q2 integration range: "
                            << "[" << Q2.min << ", " << Q2.max << "] GeV^2";
  return Q2;
}
//___________________________________________________________________________
void DISKinematicsGenerator::SetKineXY(const Interaction * interaction) const
{
// W,Q^2 => x,y using a) W^2 - M^2 = 2*Ev*M*y*(1-x) and b)  Q^2 = 2*x*y*M*Ev

  //-- get initial state information
  const InitialState & init_state = interaction->GetInitialState();
  double Ev  = init_state.GetProbeE(kRfStruckNucAtRest);
  double M   = init_state.GetTarget().StruckNucleonP4()->M();
  double M2  = TMath::Power(M,2);

  //-- get current W,Q2
  double W  = interaction->GetKinematics().W();
  double W2 = W*W;
  double Q2 = interaction->GetKinematics().Q2();

  //-- compute x,y
  double x = Q2 / (W2-M2+Q2);
  double y = (W2-M2+Q2) / (2*M*Ev);
  assert(x>0. and x<1.);
  assert(y>0. and y<1.);

  LOG("DISKinematics", pDEBUG) << "(W,Q2) => (x = "<< x<< ", y = "<< y<< ")";

  interaction->GetKinematicsPtr()->Setx(x);
  interaction->GetKinematicsPtr()->Sety(y);
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

  double max_xsec = 0.0;

  const int    NW  = 20;
  const int    NQ2 = 20;
  const double e   = 1E-6;

  LOG("DISKinematics", pDEBUG)
                << "Computing max xsec in allowed W,Q2 phase space";

  //-- Get the physical W range taking into account any user cuts
  Range1D_t W  = this->WRange(interaction);
  LOG("DISKinematics", pDEBUG)
                    << "W range = (" << W.min << ", " << W.max << ")";
  assert(W.min>0);
  const double logWmin  = TMath::Log(W.min+e);
  const double logWmax  = TMath::Log(W.max-e);
  const double dlogW    = (logWmax - logWmin)/(NW-1);

  for(int i=0; i<NW; i++) {
     double gW = TMath::Exp(logWmin + i*dlogW);
     interaction->GetKinematicsPtr()->SetW(gW);

     //-- Get the physical Q2 range (for current W) taking into account
     //   any user cuts
     Range1D_t Q2 = this->Q2Range(interaction);
     LOG("DISKinematics", pDEBUG)
        << "W = " << gW << ", Q^2 range = ("
                                 << Q2.min << ", " << Q2.max << ")";

     if(Q2.min<=0 || Q2.max <=0 || Q2.max-Q2.min<= kMinQ2Limit) continue;

     const double logQ2min = TMath::Log(Q2.min+e);
     const double logQ2max = TMath::Log(Q2.max-e);
     const double dlogQ2   = (logQ2max - logQ2min)/(NQ2-1);

     double xseclast = -1;
     bool   increasing;

     for(int j=0; j<NQ2; j++) {
         double gQ2 = TMath::Exp(logQ2min + j*dlogQ2);
         interaction->GetKinematicsPtr()->SetQ2(gQ2);

         // W,Q2 -> x,y
         this->SetKineXY(interaction);

         // update maximum xsec
         double xsec = fXSecModel->XSec(interaction);
         max_xsec = TMath::Max(xsec, max_xsec);

         increasing = xsec-xseclast>0;
         xseclast   = xsec;
         if(!increasing) break;
     } // Q2
  }// W

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  max_xsec *= fSafetyFactor;

  SLOG("DISKinematics", pDEBUG) << interaction->AsString();
  SLOG("DISKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("DISKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________

