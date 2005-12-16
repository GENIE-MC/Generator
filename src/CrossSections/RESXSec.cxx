//____________________________________________________________________________
/*!

\class    genie::RESXSec

\brief    Computed the RES Cross Section.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "CrossSections/RESXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/IntegratorI.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RESXSec::RESXSec() :
XSecAlgorithmI("genie::RESXSec")
{

}
//____________________________________________________________________________
RESXSec::RESXSec(string config) :
XSecAlgorithmI("genie::RESXSec", config)
{

}
//____________________________________________________________________________
RESXSec::~RESXSec()
{

}
//____________________________________________________________________________
double RESXSec::XSec(const Interaction * interaction) const
{
  LOG("RESXSec", pINFO) << *fConfig;

  //-- Get neutrino energy in the struck nucleon rest frame
  const InitialState & init_state = interaction -> GetInitialState();
  double Ev  = init_state.GetProbeE(kRfStruckNucAtRest);

  //-- check energy threshold
  double EvThr = utils::kinematics::EnergyThreshold(interaction);
  if(Ev <= EvThr) return 0.;

  //-- Get W integration range
  Range1D_t rW = this->WRange(interaction);

  double dW = (rW.max - rW.min) / (fNW-1);

  //-- Define the integration grid & instantiate a FunctionMap
  UnifGrid gridW;
  gridW.AddDimension(fNW, rW.min, rW.max);

  FunctionMap dxsec_dW(gridW); // dxsec/dW vs W

  //-- loop over the phase space points & compute the cross section
  for(int i = 0; i < fNW; i++) {

    double W = rW.min + i*dW;
    interaction->GetKinematicsPtr()->SetW(W);

    // Q^2 limits depend on W
    Range1D_t rQ2 = this->Q2Range(interaction);

    double d1xsec = 0.; // Q2*d^2xsec/dWdQ2 integral over logQ2

    if(TMath::Abs(rQ2.max-rQ2.min) > 1e-3) {

      double logQ2min = TMath::Log( rQ2.min );
      double logQ2max = TMath::Log( rQ2.max );
      double dlogQ2   = ( logQ2max - logQ2min) / (fNlogQ2-1);

      UnifGrid gridQ2;
      gridQ2.AddDimension(fNlogQ2, logQ2min, logQ2max);

      FunctionMap Q2d2xsec_dWdQ2(gridQ2); // Q2*d^2xsec/dWdQ2 vs Q2 [@ fixed W]

      for(int j = 0; j < fNlogQ2; j++) {

        double Q2 = TMath::Exp(logQ2min + j * dlogQ2);
        interaction->GetKinematicsPtr()->SetQ2(Q2);

        //-- compute d^2xsec / dW dQ2
        double d2xsec = fPartialXSecAlg->XSec(interaction);
        SLOG("RESXSec", pDEBUG)
              << "d^2xsec[RES]/dQ2dW (Q2 = " << Q2
                    << ", W = " << W << ", Ev = " << Ev << ") = " << d2xsec;

        //-- push Q2*(d^2xsec/dWdQ2) to the FunctionMap
        Q2d2xsec_dWdQ2.AddPoint(Q2*d2xsec, j);
      } //Q2

      //-- integrate d^2xsec/dWdQ2 over logQ2
      d1xsec = fIntegrator->Integrate(Q2d2xsec_dWdQ2);

    } // Q2min != Q2max

    SLOG("RESXSec", pDEBUG)
              << "Integral{ dlogQ2 * Q2 * d^2xsec[RES]/dQ2dW ("
                    << "W = " << W << ", Ev = " << Ev << ") } = " << d1xsec;

    //-- add integral to dxsec/dW vs W function map
    dxsec_dW.AddPoint(d1xsec, i);
  } //W

  //----- integrate dxsec/dW over W
  double xsec = fIntegrator->Integrate(dxsec_dW);

  SLOG("RESXSec", pINFO)  << "xsec[RES] (Ev = " << Ev << " GeV) = " << xsec;

  return xsec;
}
//____________________________________________________________________________
void RESXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void RESXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void RESXSec::LoadConfigData(void)
{
  //-- Get number of integration steps
  int   fNW      = fConfig->GetIntDef("nW",     31);
  int   fNlogQ2  = fConfig->GetIntDef("nLogQ2", 61);
  assert(fNW > 1 && fNlogQ2 > 1);
}
//____________________________________________________________________________
void RESXSec::LoadSubAlg(void)
{
  fPartialXSecAlg = 0;
  fIntegrator     = 0;

  //-- Get the requested d^2xsec/dxdy xsec algorithm to use
  fPartialXSecAlg =
          dynamic_cast<const XSecAlgorithmI *> (SubAlg(
                          "partial-xsec-alg-name", "partial-xsec-param-set"));

  //-- get specified integration algorithm from the config. registry
  //   for th eintegrations to follow (or use default = Simpson1D)
  string intgr = fConfig->GetStringDef("integrator", "genie::Simpson1D");

  AlgFactory * algf = AlgFactory::Instance();
  fIntegrator = dynamic_cast<const IntegratorI *> (algf->GetAlgorithm(intgr));

  assert( fPartialXSecAlg );
  assert( fIntegrator     );
}
//____________________________________________________________________________
Range1D_t RESXSec::WRange(const Interaction * interaction) const
{
  //-- Get the physically allowed W range for this interaction and allow the
  //   user inputs (if any) to narrow it

  Range1D_t rW = utils::kinematics::WRange(interaction); // physical range

  LOG("RESXSec", pDEBUG)
       << "Physical W range: " << "[" << rW.min << ", " << rW.max << "] GeV";

  // user cuts
  double Wmin = fConfig->GetDoubleDef("Wmin", -1.0);
  double Wmax = fConfig->GetDoubleDef("Wmax",  1e9);

  utils::kinematics::ApplyCutsToKineLimits(rW,  Wmin,  Wmax );

  LOG("RESXSec", pDEBUG)
       << "Physical & User W range: "
                               << "[" << rW.min << ", " << rW.max << "] GeV";
  return rW;
}
//___________________________________________________________________________
Range1D_t RESXSec::Q2Range(const Interaction * interaction) const
{
  //-- Get the physically allowed Q2 range for this interaction and allow the
  //   user inputs (if any) to narrow it

  Range1D_t rQ2 = utils::kinematics::Q2Range_W(interaction); // physical range

  LOG("RESXSec", pDEBUG) << "Physical Q2 range: "
                         << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";

  // user cuts
  double Q2min = fConfig->GetDoubleDef("Q2min", -1.0);
  double Q2max = fConfig->GetDoubleDef("Q2max",  1e9);

  utils::kinematics::ApplyCutsToKineLimits(rQ2, Q2min, Q2max );

  LOG("RESXSec", pDEBUG)
       << "Physical & User Q2 range: "
                         << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
  return rQ2;
}
//___________________________________________________________________________


