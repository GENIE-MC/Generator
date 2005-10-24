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

  //-- Get the requested d^2xsec/dxdy xsec algorithm to use

  const Algorithm * xsec_alg_base = this->SubAlg(
                           "partial-xsec-alg-name", "partial-xsec-param-set");
  const XSecAlgorithmI * partial_xsec_alg =
                         dynamic_cast<const XSecAlgorithmI *> (xsec_alg_base);

  //-- Get neutrino energy in the struck nucleon rest frame

  const InitialState & init_state = interaction -> GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double Ev  = p4->Energy();

  delete p4;

  //-- check energy threshold

  double EvThr = utils::kinematics::EnergyThreshold(interaction);
  if(Ev <= EvThr) return 0.;

  //-- Get number of integration steps

  int   nW      = this->NW();
  int   nlogQ2  = this->NLogQ2();

  assert( nW > 1 && nlogQ2 > 1 );

  //-- get specified integration algorithm from the config. registry
  //   for th eintegrations to follow (or use default = Simpson1D)

  const IntegratorI * integrator = this->Integrator();

  //-- Get W integration range

  Range1D_t rW = this->WRange(interaction);

  double dW = (rW.max - rW.min) / (nW-1);

  //-- Define the integration grid & instantiate a FunctionMap

  UnifGrid gridW;
  gridW.AddDimension(nW, rW.min, rW.max);

  FunctionMap dxsec_dW(gridW); // dxsec/dW vs W

  for(int i = 0; i < nW; i++) {

    double W = rW.min + i*dW;
    interaction->GetScatParamsPtr()->Set("W",  W );

    // Q^2 limits depend on W
    Range1D_t rQ2 = this->Q2Range(interaction);

    double d1xsec = 0.; // Q2*d^2xsec/dWdQ2 integral over logQ2

    if(TMath::Abs(rQ2.max-rQ2.min) > 1e-3) {

      double logQ2min = TMath::Log( rQ2.min );
      double logQ2max = TMath::Log( rQ2.max );
      double dlogQ2   = ( logQ2max - logQ2min) / (nlogQ2-1);

      UnifGrid gridQ2;
      gridQ2.AddDimension(nlogQ2, logQ2min, logQ2max);

      FunctionMap Q2d2xsec_dWdQ2(gridQ2); // Q2*d^2xsec/dWdQ2 vs Q2 [@ fixed W]

      for(int j = 0; j < nlogQ2; j++) {

        double Q2 = TMath::Exp(logQ2min + j * dlogQ2);
        interaction->GetScatParamsPtr()->Set("Q2", Q2);

        //-- compute d^2xsec / dW dQ2
        double d2xsec = partial_xsec_alg->XSec(interaction);

        SLOG("RESXSec", pDEBUG)
              << "d^2xsec[RES]/dQ2dW (Q2 = " << Q2
                    << ", W = " << W << ", Ev = " << Ev << ") = " << d2xsec;

        //-- push Q2*(d^2xsec/dWdQ2) to the FunctionMap

        Q2d2xsec_dWdQ2.AddPoint(Q2*d2xsec, j);
      } //Q2

      //-- integrate d^2xsec/dWdQ2 over logQ2
      d1xsec = integrator->Integrate(Q2d2xsec_dWdQ2);

    } // Q2min != Q2max

    SLOG("RESXSec", pDEBUG)
              << "Integral{ dlogQ2 * Q2 * d^2xsec[RES]/dQ2dW ("
                    << "W = " << W << ", Ev = " << Ev << ") } = " << d1xsec;

    //-- add integral to dxsec/dW vs W function map
    dxsec_dW.AddPoint(d1xsec, i);
  } //W

  //----- integrate dxsec/dW over W

  double xsec = integrator->Integrate(dxsec_dW);

  SLOG("RESXSec", pINFO)  << "xsec[RES] (Ev = " << Ev << " GeV) = " << xsec;

  return xsec;
}
//____________________________________________________________________________
const IntegratorI * RESXSec::Integrator(void) const
{
// Returns the specified (in the config. registry) integration algorithm
// If none is specified it returns a Simpson2D Integration algorithm

  string integrator_name;

  if( fConfig->Exists("integrator") )
                          fConfig->Get("integrator", integrator_name );
  else integrator_name = "genie::Simpson1D";

  //----- Get the requested algorithm from the algorithm factory

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * alg_base = algf->GetAlgorithm(integrator_name);

  const IntegratorI * integrator =
                              dynamic_cast<const IntegratorI *> (alg_base);

  assert(integrator);

  return integrator;
}
//____________________________________________________________________________
Range1D_t RESXSec::WRange(const Interaction * interaction) const
{
  //-- Get the physically allowed W range for this interaction and allow the
  //   user inputs (if any) to narrow it

  Range1D_t rW = utils::kinematics::WRange(interaction); // physical range

  LOG("RESXSec", pDEBUG)
       << "Physical W range: " << "[" << rW.min << ", " << rW.max << "] GeV";

  double Wmin = this->Wmin(); // user cuts
  double Wmax = this->Wmax();

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

  double Q2min = this->Q2min(); // user cuts
  double Q2max = this->Q2max();

  utils::kinematics::ApplyCutsToKineLimits(rQ2, Q2min, Q2max );

  LOG("RESXSec", pDEBUG)
       << "Physical & User Q2 range: "
                         << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
  return rQ2;
}
//___________________________________________________________________________
int RESXSec::NW(void) const
{
  return ( fConfig->Exists("nW") ) ? fConfig->GetInt("nW") : 31;
}
//____________________________________________________________________________
double RESXSec::Wmin(void) const
{
  return ( fConfig->Exists("Wmin") ) ? fConfig->GetDouble("Wmin") : -1.0;
}
//____________________________________________________________________________
double RESXSec::Wmax(void) const
{
  return ( fConfig->Exists("Wmax") ) ? fConfig->GetDouble("Wmax") : 1e9;
}
//____________________________________________________________________________
int RESXSec::NLogQ2(void) const
{
  return ( fConfig->Exists("nLogQ2") ) ? fConfig->GetInt("nLogQ2") : 61;
}
//____________________________________________________________________________
double RESXSec::Q2min(void) const
{
  return ( fConfig->Exists("Q2min") ) ? fConfig->GetDouble("Q2min") : -1.0;
}
//____________________________________________________________________________
double RESXSec::Q2max(void) const
{
  return ( fConfig->Exists("Q2max") ) ? fConfig->GetDouble("Q2max") : 1e9;
}
//____________________________________________________________________________


