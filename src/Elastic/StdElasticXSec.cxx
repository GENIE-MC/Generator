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
  //----- get differential cross section & integrator algorithms

  LOG("InverseMuDecay", pDEBUG) << "Getting requested diff. xsec algorithm";

  const XSecAlgorithmI * pxsec =
          dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                        "partial-xsec-alg-name", "partial-xsec-param-set"));

  LOG("InverseMuDecay", pDEBUG)
           << "Getting requested integrator algorithm (or setting default)";

  string intgr = fConfig->GetStringDef("integrator", "genie::Simpson1D");
  AlgFactory * algf = AlgFactory::Instance();
  const IntegratorI * integrator =
          dynamic_cast<const IntegratorI *> (algf->GetAlgorithm(intgr));

  //----- get initial & final state information
  const InitialState & init_state = interaction->GetInitialState();
  double E  = init_state.GetProbeE(kRfStruckNucAtRest);

  //----- integrate the differential cross section
  Range1D_t rQ2 = utils::kinematics::Q2Range_M(interaction);

  const int    nsteps  = 201;
  const double dQ2     = (rQ2.max-rQ2.min)/(nsteps-1);

  UnifGrid grid;
  grid.AddDimension(nsteps, rQ2.min, rQ2.max);

  FunctionMap fmap(grid);

  for(int i = 0; i < nsteps; i++) {

    double Q2  = rQ2.min + i * dQ2;
    interaction->GetKinematicsPtr()->SetQ2(Q2);

    double dsig_dQ2  = pxsec->XSec(interaction);
    fmap.AddPoint(dsig_dQ2, i);
  }

  double sig = integrator->Integrate(fmap);

  LOG("InverseMuDecay", pDEBUG) << "*** xsec[IMD (E=" << E << ")] = " << sig;

  return sig;
}
//____________________________________________________________________________
