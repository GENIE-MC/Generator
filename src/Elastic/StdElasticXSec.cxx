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
XSecAlgorithmI()
{
  fName = "genie::StdElasticXSec";
}
//____________________________________________________________________________
StdElasticXSec::StdElasticXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::StdElasticXSec";

  FindConfig();
}
//____________________________________________________________________________
StdElasticXSec::~StdElasticXSec()
{

}
//____________________________________________________________________________
double StdElasticXSec::XSec(const Interaction * interaction) const
{
  //----- get differential cross section & integrator algorithms

  const XSecAlgorithmI * pxsec      = this->DiffXSec();
  const IntegratorI *    integrator = this->Integrator();

  //----- get initial & final state information

  const InitialState & init_state = interaction -> GetInitialState();

  TLorentzVector * nu_p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double E  = nu_p4->Energy();

  delete nu_p4;

  //----- integrate the differential cross section

  Range1D_t rQ2 = utils::kinematics::Q2Range_M(interaction);

  const int    nsteps  = 201;
  const double dQ2     = (rQ2.max-rQ2.min)/(nsteps-1);

  UnifGrid grid;

  grid.AddDimension(nsteps, rQ2.min, rQ2.max);

  FunctionMap fmap(grid);

  for(int i = 0; i < nsteps; i++) {

    double Q2  = rQ2.min + i * dQ2;

    interaction->GetScatParamsPtr()->Set("Q2",Q2);

    double dsig_dQ2  = pxsec->XSec(interaction);

    fmap.AddPoint(dsig_dQ2, i);
  }

  double sig = integrator->Integrate(fmap);

  LOG("InverseMuDecay", pDEBUG) << "*** xsec[IMD] = " << sig;

  return sig;
}
//____________________________________________________________________________
const XSecAlgorithmI * StdElasticXSec::DiffXSec(void) const
{
 LOG("InverseMuDecay", pINFO)
                           << "Getting requested differential xsec algorithm";

 assert(fConfig->Exists("partial-xsec-alg-name") &&
                                  fConfig->Exists("partial-xsec-param-set"));

 string alg_name  = fConfig->GetString("partial-xsec-alg-name");
 string param_set = fConfig->GetString("partial-xsec-param-set");

 //-- Get the requested cross section calculator from the AlgFactory

 AlgFactory * algf = AlgFactory::Instance();

 const Algorithm * alg_base = algf->GetAlgorithm(alg_name, param_set);

 const XSecAlgorithmI * xs = dynamic_cast<const XSecAlgorithmI *> (alg_base);

 assert(xs);

 return xs;
}
//____________________________________________________________________________
const IntegratorI * StdElasticXSec::Integrator(void) const
{
 LOG("InverseMuDecay", pINFO) << "Getting requested integrator algorithm";

 // read integrator name from config or set default
 string integrator_name = ( fConfig->Exists("integrator-name") ) ?
                  fConfig->GetString("integrator-name") : "genie::Simpson1D";

 // get integrator algorithm

 AlgFactory * algf = AlgFactory::Instance();

 const Algorithm * alg_base = algf->GetAlgorithm(integrator_name);

 const IntegratorI * intgr = dynamic_cast<const IntegratorI *> (alg_base);

 assert(intgr);

 return intgr;
}
//____________________________________________________________________________
