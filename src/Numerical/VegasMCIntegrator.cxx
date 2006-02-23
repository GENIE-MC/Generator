//____________________________________________________________________________
/*!

\class    genie::VegasMCIntegrator

\brief    The VEGAS (P.Lepage) adaptive MC integration algorithm based on
          importance and stratified sampling.

\ref      Numerical Recipes in C, Cambridge Univ. Press, 2002, page 316-328

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  February 23, 2006

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Numerical/VegasMCIntegrator.h"
#include "Numerical/GSFunc.h"
#include "Numerical/RandomGen.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
VegasMCIntegrator::VegasMCIntegrator():
IntegratorI("genie::VegasMCIntegrator")
{

}
//____________________________________________________________________________
VegasMCIntegrator::VegasMCIntegrator(string config) :
IntegratorI("genie::VegasMCIntegrator", config)
{

}
//____________________________________________________________________________
VegasMCIntegrator::~VegasMCIntegrator()
{

}
//____________________________________________________________________________
double VegasMCIntegrator::Integrate(GSFunc & gsfunc) const
{
  unsigned int ndim = gsfunc.NParams();

  LOG("VegasMCIntegrator", pINFO)
                   << "VEGAS MC integration in: "  << ndim << " dimensions";

  return 0;
}
//____________________________________________________________________________
void VegasMCIntegrator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  //...
  //...
}
//____________________________________________________________________________
void VegasMCIntegrator::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  //...
  //...
}
//____________________________________________________________________________

