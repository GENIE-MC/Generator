//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the new Coherent package from its previous location  (Elastic
   package)

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Physics/Coherent/XSection/DeVriesFormFactor.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;
using namespace genie::utils;

//____________________________________________________________________________
DeVriesFormFactor::DeVriesFormFactor() :
XSecAlgorithmI("genie::DeVriesFormFactor")
{

}
//____________________________________________________________________________
DeVriesFormFactor::DeVriesFormFactor(string config) :
XSecAlgorithmI("genie::DeVriesFormFactor", config)
{

}
//____________________________________________________________________________
DeVriesFormFactor::~DeVriesFormFactor()
{

}
//____________________________________________________________________________
double DeVriesFormFactor::ProtonFF( double Q ) const {

  
}
//____________________________________________________________________________
double DeVriesFormFactor::NeutronFF( double Q ) const ;
//____________________________________________________________________________
void DeVriesFormFactor::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeVriesFormFactor::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeVriesFormFactor::LoadConfig(void)
{


}
//____________________________________________________________________________
