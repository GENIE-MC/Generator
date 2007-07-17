//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - February 23, 2006

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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

