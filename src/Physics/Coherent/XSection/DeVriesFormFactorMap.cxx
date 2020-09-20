//____________________________________________________________________________
/*
  Copyright (c) 2003-2019, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Author: Marco Roda
  University of Liverpool
  <mroda \at liverpool.ac.uk>

  For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <vector>

#include "Framework/Messenger/Messenger.h"

#include "Physics/Coherent/XSection/DeVriesFormFactorMap.h"
#include "Framework/Registry/RegistryItemTypeDef.h"

using namespace genie;


DeVriesFormFactorMap::DeVriesFormFactorMap() :
  COHFormFactorI("genie::DeVriesFormFactorMap")
{

}
//____________________________________________________________________________
DeVriesFormFactorMap::DeVriesFormFactorMap(string config) :
  COHFormFactorI("genie::DeVriesFormFactorMap", config)
{

}
//____________________________________________________________________________
DeVriesFormFactorMap::DeVriesFormFactorMap( string name, string config ) :
  COHFormFactorI(name, config)
{

}
//____________________________________________________________________________
DeVriesFormFactorMap::~DeVriesFormFactorMap()
{

}
//____________________________________________________________________________
double DeVriesFormFactorMap::ProtonFF( double Q, int pdg ) const {

  const std::map<int, const genie::DeVriesFormFactor *>::const_iterator it =
    fNuclearFFs.find( pdg ) ;

  if ( it == fNuclearFFs.end() ) return 0. ;

  return it -> second -> Calculator().FormFactor( Q ) ;

}
//____________________________________________________________________________
bool DeVriesFormFactorMap::HasNucleus( int pdg ) const {

  return (fNuclearFFs.count( pdg ) > 0) ;
}
//____________________________________________________________________________
void DeVriesFormFactorMap::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeVriesFormFactorMap::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeVriesFormFactorMap::LoadConfig(void)
{

  fNuclearFFs.clear() ;

  bool good_configuration = true ;

  // read the vector of algos from the xml file
  std::vector<RgKey> keys ;
  GetParamVectKeys( "COH-DV-FormFactor", keys ) ;

  // Store pointers to subalgos in the local map
  for ( unsigned int i = 0 ; i < keys.size() ; ++i ) {

    const Algorithm * algo = SubAlg( keys[i] ) ;

    const DeVriesFormFactor * ff = dynamic_cast< const DeVriesFormFactor * >( algo ) ;

    if ( ! ff ) {
      good_configuration = false ;
      LOG("DeVriesFormFactorMap", pERROR ) << "SubAlgo with key " << keys[i] << " not retrieved" ;

    }

    if ( fNuclearFFs.count( ff -> NucleusPDG() ) > 0 ) {

      good_configuration = false ;
      LOG("DeVriesFormFactorMap", pERROR ) << "Attempt to add a second DeVries form factor for PDG " << ff -> NucleusPDG() ;
    }

    fNuclearFFs[ ff -> NucleusPDG() ] = ff ;

  }  // loop over subalgo

  if ( ! good_configuration ) {
    LOG("DeVriesFormFactorMap", pFATAL ) << "Configuration not good, exiting" ;
    exit ( 78 ) ;
  }

}
//____________________________________________________________________________
