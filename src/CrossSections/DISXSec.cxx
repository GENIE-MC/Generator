//____________________________________________________________________________
/*!

\class    genie::DISXSec

\brief    Computes the DIS Cross Section.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "CrossSections/DISXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/IntegratorI.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISXSec::DISXSec() :
XSecAlgorithmI("genie::DISXSec")
{

}
//____________________________________________________________________________
DISXSec::DISXSec(string config) :
XSecAlgorithmI("genie::DISXSec", config)
{

}
//____________________________________________________________________________
DISXSec::~DISXSec()
{

}
//____________________________________________________________________________
double DISXSec::XSec(const Interaction * interaction) const
{
  //-- Get neutrino energy in the struck nucleon rest frame
  const InitialState & init_state = interaction -> GetInitialState();
  double Ev = init_state.GetProbeE(kRfStruckNucAtRest);

  //-- Check the energy threshold
  double Ethr = utils::kinematics::EnergyThreshold(interaction);
  if(Ev <= Ethr) {
     LOG("DISXSec", pINFO) << "E = " << Ev << " <= Ethreshold = " << Ethr;
     return 0.;
  }

  //-- Define the integration grid & instantiate a FunctionMap
  UnifGrid grid;
  grid.AddDimension(fNlogx, fLogXmin, fLogXmax);
  grid.AddDimension(fNlogy, fLogYmin, fLogYmax);

  FunctionMap xyd2xsec(grid);

  //----- Loop over x,y & compute the differential xsec
  LOG("DISXSec", pDEBUG)
       << "integration limits: x = (" << fXmin << ", " << fXmax << ") "
                           << "y = (" << fYmin << ", " << fYmax << ")";

  for(int ix = 0; ix < fNlogx; ix++) {
    double x  = TMath::Exp(fLogXmin + ix * fdLogX);

    for(int iy = 0; iy < fNlogy; iy++) {
       double y = TMath::Exp(fLogYmin + iy * fdLogY);

       //-- update the scattering parameters
       interaction->GetKinematicsPtr()->Setx(x);
       interaction->GetKinematicsPtr()->Sety(y);
       bool is_physical = this->IsWithinIntegrationRange(interaction);

       double pxsec = 0.;
       if (is_physical) {
          //-- compute d^2xsec/dxdy
          pxsec = fPartialXSecAlg->XSec(interaction);
          LOG("DISXSec", pDEBUG)
              << "dxsec/dxdy (x = " << x << ", y = " << y
                          << ", Ev = " << Ev << ") = " << pxsec;
       }
       //-- push x*y*(d^2xsec/dxdy) to the FunctionMap
       xyd2xsec.AddPoint(x*y*pxsec, ix, iy);
    } //y
  } //x

  //----- Perform the numerical integration
  LOG("DISXSec", pDEBUG) << "Performing numerical integration";
  double xsec = fIntegrator->Integrate(xyd2xsec);

  LOG("DISXSec", pDEBUG)  << "xsec_dis (E = " << Ev << " GeV) = " << xsec;

  return xsec;
}
//____________________________________________________________________________
void DISXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void DISXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void DISXSec::LoadConfigData(void)
{
  //-- Get x,y from config (if they exist) or set defaults
  fNlogx = fConfig -> GetIntDef    ("n-log-x", 61   );
  fNlogy = fConfig -> GetIntDef    ("n-log-y", 61   );
  fXmin  = fConfig -> GetDoubleDef ("x-min",   0.001);
  fXmax  = fConfig -> GetDoubleDef ("x-max",   0.999);
  fYmin  = fConfig -> GetDoubleDef ("y-min",   0.001);
  fYmax  = fConfig -> GetDoubleDef ("y-max",   0.999);

  //-- Check that x,y range is meaningful
  assert(fNlogx > 1 && fNlogy > 1);
  assert(fXmax > fXmin && fXmax < 1 && fXmin < 1 && fXmax > 0 & fXmin > 0);
  assert(fYmax > fYmin && fYmax < 1 && fYmin < 1 && fYmax > 0 & fYmin > 0);

  //-- Define the integration area
  fLogXmax = TMath::Log(fXmax);
  fLogXmin = TMath::Log(fXmin);
  fLogYmax = TMath::Log(fYmax);
  fLogYmin = TMath::Log(fYmin);
  fdLogX   = (fLogXmax - fLogXmin) / (fNlogx-1);
  fdLogY   = (fLogYmax - fLogYmin) / (fNlogy-1);

  // check whether the user imposes kinematic cuts
  fWmin  = fConfig->GetDoubleDef( "Wmin",  -1  );
  fWmax  = fConfig->GetDoubleDef( "Wmax",  1e9 );
  fQ2min = fConfig->GetDoubleDef( "Q2min", -1  );
  fQ2max = fConfig->GetDoubleDef( "Q2max", 1e9 );
}
//____________________________________________________________________________
void DISXSec::LoadSubAlg(void)
{
  fPartialXSecAlg = 0;
  fIntegrator     = 0;

  //-- Get the requested d^2xsec/dxdy xsec algorithm to use
  fPartialXSecAlg =
         dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fPartialXSecAlg);

  LOG("DISXSec", pDEBUG) << *fPartialXSecAlg;

  //-- get specified integration algorithm from the config. registry
  //   or use Simpson2D if none is defined
  string intgr = fConfig->GetStringDef("integrator", "genie::Simpson2D");

  AlgFactory * algf = AlgFactory::Instance();
  fIntegrator = dynamic_cast<const IntegratorI *> (algf->GetAlgorithm(intgr));
  assert(fIntegrator);
}
//____________________________________________________________________________
bool DISXSec::IsWithinIntegrationRange(const Interaction * interaction) const
{
  //-- Allows the user to set integration limits on W, Q2. For each x,y
  //   pair, the corresponding W,Q2 pair is calculated. The physical
  //   W, Q2 range is calculated and the user-cuts taken into account for
  //   narrowing it down.

  // get physical integration range for W and Q2
  Range1D_t rW  = utils::kinematics::WRange     (interaction);
  Range1D_t rQ2 = utils::kinematics::Q2Range_xy (interaction);

  LOG("DISXSec", pDEBUG)
       << "\n Physical integration range: "
             << "W = [" << rW.min << ", " << rW.max << "] GeV "
                      << "Q2 = [" << rQ2.min << ", " << rQ2.max << "] GeV^2";

  // apply cuts
  utils::kinematics::ApplyCutsToKineLimits(rW,  fWmin,  fWmax );
  utils::kinematics::ApplyCutsToKineLimits(rQ2, fQ2min, fQ2max);

  LOG("DISXSec", pDEBUG)
       << "\n Physical && User integration range: "
            << "W = [" << rW.min << ", " << rW.max << "] GeV "
                 << "Q2 = [" << rQ2.min << ", " << rQ2.max << "] GeV^2";
  // current W, Q2

  const InitialState & init_state = interaction -> GetInitialState();

  double x  = interaction->GetKinematics().x();
  double y  = interaction->GetKinematics().y();
  double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
  double M  = init_state.GetTarget().StruckNucleonMass();
  double M2 = M*M;

  double currW2 = TMath::Max(0., M2 + 2*Ev*M*y*(1-x));
  double currW  = TMath::Sqrt(currW2);
  double currQ2 = TMath::Max(0., 2*M*Ev*x*y);

  bool in_range = utils::math::IsWithinLimits(currQ2, rQ2)
                                && utils::math::IsWithinLimits(currW, rW);

  if(!in_range) {
      LOG("DISXSec", pDEBUG) << "*** excluding from phase space: ";
  } else {
      LOG("DISXSec", pDEBUG) << "*** including in phase space: ";
  }
  LOG("DISXSec", pDEBUG)
          << "x = " << x << ", y = " << y << ", Q2 = " << currQ2
                                << ", W = " << currW << ", Ev = " << Ev;
  return in_range;
}
//____________________________________________________________________________
