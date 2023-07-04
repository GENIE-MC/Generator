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

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"

#include "Framework/Utils/StringUtils.h"

using namespace genie;


//____________________________________________________________________________
DeVriesFormFactor::DeVriesFormFactor() :
  Algorithm("genie::DeVriesFormFactor"),
  fCalculator( {}, 0., 0., 0. ),
  fPDG(0) { ; }
//____________________________________________________________________________
DeVriesFormFactor::DeVriesFormFactor(string config) :
  Algorithm("genie::DeVriesFormFactor", config),
  fCalculator( {}, 0., 0., 0. ),
  fPDG(0) { ; }
//____________________________________________________________________________
DeVriesFormFactor::~DeVriesFormFactor()
{

}
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

  // load coefficients
  vector<double> cs;
  this -> GetParamVect( "DV-Coefficient", cs ) ;

  double r;
  GetParam( "DV-Radius", r ) ;
  r *= units::fm ;

  double Qmax, Qmin;
  GetParamDef( "DV-QMin", Qmin, 0. ) ;
  GetParam( "DV-QMax", Qmax ) ;
  // Qmax and Qmin are released in the DeVries paper in fm. 
  // and so is the input in the xml file
  // So we need to convert it back to natural units
  Qmin /= units::fm ;
  Qmax /= units::fm ;
  
  fCalculator = FourierBesselFFCalculator( cs, r, Qmin, Qmax ) ;

  GetParam( "DV-Nucleus", fPDG ) ;

  LOG("DeVriesFormFactor", pINFO) << "Loaded " << cs.size() << " coefficients for nucleus " << fPDG ;

}
//____________________________________________________________________________
