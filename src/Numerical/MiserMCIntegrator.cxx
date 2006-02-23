//____________________________________________________________________________
/*!

\class    genie::MiserMCIntegrator

\brief    The MISER adaptive MC integration algorithm based on recursive
          stratified sampling.

\ref      Numerical Recipes in C, Cambridge Univ. Press, 2002, page 316-328

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  February 23, 2006

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Numerical/MiserMCIntegrator.h"
#include "Numerical/GSFunc.h"
#include "Numerical/RandomGen.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
MiserMCIntegrator::MiserMCIntegrator():
IntegratorI("genie::MiserMCIntegrator")
{

}
//____________________________________________________________________________
MiserMCIntegrator::MiserMCIntegrator(string config) :
IntegratorI("genie::MiserMCIntegrator", config)
{

}
//____________________________________________________________________________
MiserMCIntegrator::~MiserMCIntegrator()
{

}
//____________________________________________________________________________
double MiserMCIntegrator::Integrate(GSFunc & gsfunc) const
{
  unsigned int ndim = gsfunc.NParams();

  LOG("MiserMCIntegrator", pINFO)
                   << "VEGAS MC integration in: "  << ndim << " dimensions";

  return 0;
}
//____________________________________________________________________________
void MiserMCIntegrator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  //...
  //...
}
//____________________________________________________________________________
void MiserMCIntegrator::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  //...
  //...
}
//____________________________________________________________________________

