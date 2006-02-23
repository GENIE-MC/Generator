//____________________________________________________________________________
/*!

\class    genie::StdElasticXSec

\brief    Standard v+N / vbar+N elastic scattering cross section.

          StdElasticPXSec is a concrete implementation of the
          XSecAlgorithmI interface. \n

\ref      L.A.Ahrens et al., Physical Review D, VOL 35,3:785 (1987)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 15, 2005

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/GXSecFunc.h"
#include "Elastic/StdElasticXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
StdElasticXSec::StdElasticXSec() :
XSecAlgorithmI("genie::StdElasticXSec")
{

}
//____________________________________________________________________________
StdElasticXSec::StdElasticXSec(string config) :
XSecAlgorithmI("genie::StdElasticXSec", config)
{

}
//____________________________________________________________________________
StdElasticXSec::~StdElasticXSec()
{

}
//____________________________________________________________________________
double StdElasticXSec::XSec(const Interaction * interaction) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get initial & final state information
  const InitialState & init_state = interaction->GetInitialState();
  double E  = init_state.GetProbeE(kRfStruckNucAtRest);

  //----- integrate the differential cross section
  Range1D_t rQ2 = utils::kinematics::Q2Range_M(interaction);

  GXSecFunc * func = new Integrand_DXSec_DQ2_E(fDiffXSecModel, interaction);
  func->SetParam(0,"Q2",rQ2);
  double xsec = fIntegrator->Integrate(*func);

  LOG("Elastic", pDEBUG) << "*** XSec[EL] (E=" << E << ") = " << xsec;

  delete func;
  return xsec;
}
//____________________________________________________________________________
bool StdElasticXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  return true;
}
//____________________________________________________________________________
bool StdElasticXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  return true;
}
//____________________________________________________________________________
void StdElasticXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void StdElasticXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void StdElasticXSec::LoadConfig(void)
{
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed and sets configuration variables to avoid looking up at the Registry
// all the time.

  //-- get an algorithm to calculate differential cross sections dxsec/dQ2
  fDiffXSecModel =
          dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fDiffXSecModel);

  //-- get the specified integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *> (
                 this->SubAlg("integrator-alg-name", "integrator-param-set"));
  assert(fIntegrator);
}
//____________________________________________________________________________

