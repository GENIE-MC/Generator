//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - February 23, 2006

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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

