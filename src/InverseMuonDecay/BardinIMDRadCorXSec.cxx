//____________________________________________________________________________
/*!

\class    genie::BardinIMDRadCorXSec

\brief    Computes the Inverse Muon Decay cross section using the Bardin -
          Dokuchaeva model which includes all 1-loop radiative corrections. \n

          This algorithm merely integrates the Bardin differential IMD cross
          section. The specific differential cross section algorithm is
          specified in this algorithm's XML config file.

          The exact 'type' of the cross section depends on the specified
          differential IMD cross section algorithm. It can be a 'trully'
          inclusive IMD cross section or a cross section where part of the
          brem cross section contribution is subtracted
          (for futher details, see the documentation of the Bardin-Dokuchaeva
          model diffential IMD cross section algorithms and the Bardin paper,
          cited below).

          BardinIMDRadCorXSec is a concrete implementation of the
          XSecAlgorithmI interface. \n

\ref      D.Yu.Bardin and V.A.Dokuchaeva, Nucl.Phys.B287:839 (1987)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 14, 2005

*/
//____________________________________________________________________________

#include <iostream>

#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "InverseMuonDecay/BardinIMDRadCorXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BardinIMDRadCorXSec::BardinIMDRadCorXSec() :
XSecAlgorithmI("genie::BardinIMDRadCorXSec")
{

}
//____________________________________________________________________________
BardinIMDRadCorXSec::BardinIMDRadCorXSec(string config) :
XSecAlgorithmI("genie::BardinIMDRadCorXSec", config)
{

}
//____________________________________________________________________________
BardinIMDRadCorXSec::~BardinIMDRadCorXSec()
{

}
//____________________________________________________________________________
double BardinIMDRadCorXSec::XSec(const Interaction * interaction) const
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
  LOG("InverseMuDecay", pDEBUG) << "*** xsec[IMD] = " << xsec;
  return xsec;
}
//____________________________________________________________________________
bool BardinIMDRadCorXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  return true;
}
//____________________________________________________________________________
bool BardinIMDRadCorXSec::ValidKinematics(
                                        const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction -> GetInitialState();

  double E = init_state.GetProbeE(kRfLab);
  double s = kElectronMass_2 + 2*kElectronMass*E;

  //-- check if it is kinematically allowed
  if(s < kMuonMass_2) {
     LOG("InverseMuDecay", pINFO)
        << "Ev = " << E << " (s = " << s << ") is below threshold (s-min = "
        << kMuonMass_2 << ") for IMD";
     return false;
  }
  return true;
}
//____________________________________________________________________________
void BardinIMDRadCorXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BardinIMDRadCorXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BardinIMDRadCorXSec::LoadConfig(void)
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
  LOG("QELXSec", pDEBUG) << "Number of integration (logQ2) bins = " << fNBins;
  assert(fNBins>2);
}
//____________________________________________________________________________


