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
  //-- get a differential IMD cross section sub-algorithm
  const XSecAlgorithmI * pxsec =
           dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                          "partial-xsec-alg-name", "partial-xsec-param-set"));

  //-- get a 1-D nuerical integrator or set default
  const IntegratorI * integrator =
          dynamic_cast<const IntegratorI *> (this->SubAlgWithDefault(
                                 "integrator-name","","genie::Simpson1D",""));
  const int    nsteps  = 201;
  const double min     = 0;
  const double max     = 0.999;
  const double step    = (max-min)/(nsteps-1);

  UnifGrid grid;
  grid.AddDimension(nsteps, min, max);

  FunctionMap fmap(grid);

  //-- all kinematical cuts (energy threshold, physical y range) are
  //   applied within the differential cross section algorithm - returns 0
  //   if kinematic params are not valid.

  for(int i = 0; i < nsteps; i++) {
    double y = min + i * step;
    interaction->GetKinematicsPtr()->Sety(y);
    double dsig_dy = pxsec->XSec(interaction);
    fmap.AddPoint(dsig_dy, i);
  }
  double sig = integrator->Integrate(fmap);

  LOG("InverseMuDecay", pDEBUG) << "*** xsec[IMD] = " << sig;

  return sig;
}
//____________________________________________________________________________
