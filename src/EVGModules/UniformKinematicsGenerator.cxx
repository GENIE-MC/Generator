//____________________________________________________________________________
/*!

\class   genie::UniformKinematicsGenerator

\brief   Generates values for the kinematic variables describing QEL,RES,DIS
         or COH neutrino interaction events (depending on the interaction
         already found at the event record it visits) using a flat probability
         distribution.

         Only use this one if you want to generate samples with many events
         from improbable phase space regions without having to generate very
         large statistics.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 29, 2005

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "EVGModules/UniformKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//___________________________________________________________________________
UniformKinematicsGenerator::UniformKinematicsGenerator() :
EventRecordVisitorI("genie::UniformKinematicsGenerator")
{

}
//___________________________________________________________________________
UniformKinematicsGenerator::UniformKinematicsGenerator(string config) :
EventRecordVisitorI("genie::UniformKinematicsGenerator", config)
{

}
//___________________________________________________________________________
UniformKinematicsGenerator::~UniformKinematicsGenerator()
{

}
//___________________________________________________________________________
void UniformKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Selects kinematic variables and adds them to the event record's summary.
// This is a 'Uniform Kinematics Generator', so no differential xsec model
// is attached. Every point in the available phase space is equally probable.
// It is a generic visitor that would work for any interaction.

  Interaction * interaction = evrec->GetInteraction();
  const ProcessInfo & proc_info = interaction->GetProcessInfo();

  if(proc_info.IsQuasiElastic())
         this->GenerateUnifQELKinematics(evrec);
  else if (proc_info.IsResonant())
         this->GenerateUnifRESKinematics(evrec);
  else if (proc_info.IsDeepInelastic())
         this->GenerateUnifDISKinematics(evrec);
  else if (proc_info.IsCoherent())
         this->GenerateUnifCOHKinematics(evrec);
  else {
    LOG("UnifKinematics", pERROR)
       << "Can not generate uniform kinematics for scattering type:"
                                     << proc_info.ScatteringTypeAsString();
  }
}
//___________________________________________________________________________
void UniformKinematicsGenerator::GenerateUnifQELKinematics(
                                                    GHepRecord * evrec) const
{
  Interaction * interaction = evrec->GetInteraction();

  // compute physical Q^2 range
  Range1D_t Q2 = utils::kinematics::Q2Range_M(interaction);

  // generate/set a Q^2 in available phase space (const probability)
  RandomGen * rnd = RandomGen::Instance();
  double gQ2 = Q2.min + (Q2.max-Q2.min) * rnd->Random1().Rndm();
  interaction->GetKinematicsPtr()->SetQ2(gQ2);

  LOG("UnifKinematics", pINFO) << "Selected: Q^2 = " << gQ2;

  evrec->SetDiffXSec(0);
}
//___________________________________________________________________________
void UniformKinematicsGenerator::GenerateUnifRESKinematics(
                                                    GHepRecord * evrec) const
{
  Interaction * interaction = evrec->GetInteraction();

  // compute physical W range
  Range1D_t W  = utils::kinematics::WRange(interaction);

  bool found = false;
  register unsigned int iter = 0;
  double gW = -1, gQ2 = -1;

  RandomGen * rnd = RandomGen::Instance();

  while(!found) {

     // generate/set an invariant mass (const probability)
     gW = W.min + (W.max - W.min) * rnd->Random1().Rndm();
     interaction->GetKinematicsPtr()->SetW(gW);

     // compute physical Q^2 for selected W
     Range1D_t Q2 = utils::kinematics::Q2Range_W(interaction);

     if(Q2.min < Q2.max) {
        // generate/set a Q^2 (const probability)
        gQ2 = Q2.min + (Q2.max-Q2.min) * rnd->Random1().Rndm();
        interaction->GetKinematicsPtr()->SetQ2(gQ2);
        found = true;
     }

     iter++;
     if(iter > kRjMaxIterations) {
         LOG("UnifKinematics", pERROR)
             << "*** Could not select a valid (W,Q^2) pair after "
                                                  << iter << " iterations";
         exit(1);
     }
  }

  LOG("UnifKinematics", pINFO) << "Selected: W   = " << gW;
  LOG("UnifKinematics", pINFO) << "Selected: Q^2 = " << gQ2;

  evrec->SetDiffXSec(0);
}
//___________________________________________________________________________
void UniformKinematicsGenerator::GenerateUnifDISKinematics(
                                                    GHepRecord * evrec) const
{
  Interaction * interaction = evrec->GetInteraction();

  // get initial state information
  const InitialState & init_state = interaction->GetInitialState();

  double Ev  = init_state.GetProbeE(kRfStruckNucAtRest);
  double M   = init_state.GetTarget().StruckNucleonP4()->M(); // can be off m/shell
  double M2  = TMath::Power(M,2);

  // select kinematics

  bool found = false;
  register unsigned int iter = 0;

  double gx = -1, gy = -1, gW = -1, gQ2 = -1;

  RandomGen * rnd = RandomGen::Instance();

  while(!found) {
     // generate x,y in [0,1]
     gx = rnd->Random1().Rndm();
     gy = rnd->Random1().Rndm();

     interaction->GetKinematicsPtr()->Setx(gx);
     interaction->GetKinematicsPtr()->Sety(gy);

     // check whether generated x,y correspond to valid kinematics
     gW  = TMath::Sqrt( utils::math::NonNegative(M2+2*Ev*M*gy*(1-gx)) );
     gQ2 = 2*gx*gy*M*Ev;

     //-- get the physical W, Q^2 range
     Range1D_t W  = utils::kinematics::WRange(interaction);
     Range1D_t Q2 = utils::kinematics::Q2Range_xy(interaction);

     bool valid_W  = utils::math::IsWithinLimits( gW, W  );
     bool valid_Q2 = utils::math::IsWithinLimits( gQ2,Q2 );

     found = (valid_W && valid_Q2);

     iter++;
     if(iter > kRjMaxIterations) {
         LOG("UnifKinematics", pERROR)
             << "*** Could not select a valid (W,Q^2) pair after "
                                                  << iter << " iterations";
         exit(1);
     }
  }

  LOG("UnifKinematics", pINFO) << "Selected: x   = " << gx;
  LOG("UnifKinematics", pINFO) << "Selected: y   = " << gy;
  LOG("UnifKinematics", pINFO) << "Selected: W   = " << gW;
  LOG("UnifKinematics", pINFO) << "Selected: Q^2 = " << gQ2;

  evrec->SetDiffXSec(0);
}
//___________________________________________________________________________
void UniformKinematicsGenerator::GenerateUnifCOHKinematics(
                                                   GHepRecord * evrec) const
{
  LOG("UnifKinematics", pWARN) << "Kinematics generator not implemented!!";

  Interaction * interaction = evrec->GetInteraction();

  double Ev  = interaction->GetInitialState().GetProbeE(kRfLab);
  double Mpi = kPionMass;

  double ymin = Mpi/Ev;
  double ymax = 1.;

  RandomGen * rnd = RandomGen::Instance();
  double gx = rnd->Random1().Rndm(); // x in [0,1]
  double gy = ymin + (ymax-ymin)*rnd->Random1().Rndm(); // y in [ymin, ymax]

  LOG("UnifKinematics", pINFO) << "Selected: x = " << gx;
  LOG("UnifKinematics", pINFO) << "Selected: y = " << gy;

  interaction->GetKinematicsPtr()->Setx(gx);
  interaction->GetKinematicsPtr()->Sety(gy);
  evrec->SetDiffXSec(0);
}
//___________________________________________________________________________

