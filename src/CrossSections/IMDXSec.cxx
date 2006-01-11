//____________________________________________________________________________
/*!

\class    genie::IMDXSec

\brief    Computes the Inverse Muon Decay cross section

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 14, 2005

*/
//____________________________________________________________________________

#include <iostream>

#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/IMDXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
IMDXSec::IMDXSec() :
XSecAlgorithmI("genie::IMDXSec")
{

}
//____________________________________________________________________________
IMDXSec::IMDXSec(string config) :
XSecAlgorithmI("genie::IMDXSec", config)
{

}
//____________________________________________________________________________
IMDXSec::~IMDXSec()
{

}
//____________________________________________________________________________
double IMDXSec::XSec(const Interaction * interaction) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const double e       = 1e-4;
  const double ymin    = e;
  const double ymax    = 1-e;
  const double logymax = TMath::Log(ymax);
  const double logymin = TMath::Log(ymin);
  const double dlogy   = (logymax-logymin)/(fNBins-1);

  UnifGrid grid;
  grid.AddDimension(fNBins, ymin, ymax);

  FunctionMap fmap(grid);

  //-- all kinematical cuts (energy threshold, physical y range) are
  //   applied within the differential cross section algorithm - returns 0
  //   if kinematic params are not valid.

  for(int i = 0; i < fNBins; i++) {
    double y = TMath::Exp(logymin + i * dlogy);
    interaction->GetKinematicsPtr()->Sety(y);

    double dsig_dy = fDiffXSecModel->XSec(interaction);
    fmap.AddPoint(dsig_dy, i);
  }

  double xsec = fIntegrator->Integrate(fmap);
  LOG("IMD", pDEBUG) << "*** xsec[IMD] = " << xsec;
  return xsec;
}
//____________________________________________________________________________
bool IMDXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  return true;
}
//____________________________________________________________________________
bool IMDXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction -> GetInitialState();

  double E = init_state.GetProbeE(kRfLab);
  double s = kElectronMass_2 + 2*kElectronMass*E;

  //-- check if it is kinematically allowed
  if(s < kMuonMass_2) {
     LOG("IMD", pINFO)
        << "Ev = " << E << " (s = " << s << ") is below threshold (s-min = "
        << kMuonMass_2 << ") for IMD";
     return false;
  }
  return true;
}
//____________________________________________________________________________
void IMDXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void IMDXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void IMDXSec::LoadConfig(void)
{
  fDiffXSecModel = 0;
  fIntegrator    = 0;

  //----- Get an algorithm to calculate differential cross sections dxsec/dQ2
  fDiffXSecModel =
          dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fDiffXSecModel);

  //----- Get an integrator
  fIntegrator =
          dynamic_cast<const IntegratorI *> (this->SubAlgWithDefault(
                                 "integrator-name","","genie::Simpson1D",""));
  assert(fIntegrator);

  //----- get the input number of dxsec/dlogQ2 points for num. integration
  //      or use a default if no number is specified
  //      (must be odd number for the Simpson rule)
  fNBins = fConfig->GetIntDef("N-logQ2-bins", 101);
  LOG("IMD", pDEBUG) << "Number of integration (logQ2) bins = " << fNBins;
  assert(fNBins>2);
}
//____________________________________________________________________________


