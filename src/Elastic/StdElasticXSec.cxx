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

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Elastic/StdElasticXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
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

  const double logQ2min = TMath::Log(rQ2.min);
  const double logQ2max = TMath::Log(rQ2.max);
  const double dlogQ2   = (logQ2max-logQ2min)/(fNBins-1);

  UnifGrid grid;
  grid.AddDimension(fNBins, logQ2min, logQ2max);

  FunctionMap fmap(grid);

  for(int i = 0; i < fNBins; i++) {

    double Q2  = TMath::Exp(logQ2min + i * dlogQ2);
    interaction->GetKinematicsPtr()->SetQ2(Q2);

    double dsig_dQ2  = fDiffXSecModel->XSec(interaction);
    fmap.AddPoint(Q2*dsig_dQ2, i);
  }

  double xsec = fIntegrator->Integrate(fmap);
  LOG("Elastic", pDEBUG) << "*** XSec[EL] (E=" << E << ") = " << xsec;
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

  //----- Get an algorithm to calculate differential cross sections dxsec/dQ2
  fDiffXSecModel =
          dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fDiffXSecModel);

  //----- Get an integrator
  string intgr = fConfig->GetStringDef("integrator", "genie::Simpson1D");

  AlgFactory * algf = AlgFactory::Instance();
  fIntegrator = dynamic_cast<const IntegratorI *>(algf->GetAlgorithm(intgr));

  assert(fIntegrator);

  //----- get the input number of dxsec/dlogQ2 points for num. integration
  //      or use a default if no number is specified
  //      (must be odd number for the Simpson rule)
  fNBins = fConfig->GetIntDef("N-logQ2-bins", 131);
  LOG("Elastic", pDEBUG) << "Number of integration (logQ2) bins = " << fNBins;
  assert(fNBins>2);
}
//____________________________________________________________________________

