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

#include "Physics/Coherent/XSection/COHFormFactorMap.h"
#include "Framework/Registry/RegistryItemTypeDef.h"

using namespace genie;


COHFormFactorMap::COHFormFactorMap() :
  COHFormFactorI("genie::COHFormFactorMap")
{

}
//____________________________________________________________________________
COHFormFactorMap::COHFormFactorMap(string config) :
  COHFormFactorI("genie::COHFormFactorMap", config)
{

}
//____________________________________________________________________________
COHFormFactorMap::COHFormFactorMap( string name, string config ) :
  COHFormFactorI(name, config)
{

}
//____________________________________________________________________________
COHFormFactorMap::~COHFormFactorMap()
{

}
//____________________________________________________________________________
double COHFormFactorMap::ProtonFF( double Q, int pdg ) const {

  const std::map<int, const genie::DeVriesFormFactor *>::const_iterator it =
    fNuclearFFs.find( pdg ) ;

  if ( it == fNuclearFFs.end() ) return 0. ;

  return it -> second -> Calculator().FormFactor( Q ) ;

}
//____________________________________________________________________________
bool COHFormFactorMap::HasNucleus( int pdg ) const {

  return (fNuclearFFs.count( pdg ) > 0) ;
}
//____________________________________________________________________________
void COHFormFactorMap::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHFormFactorMap::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHFormFactorMap::LoadConfig(void)
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
      LOG("COHFormFactorMap", pERROR ) << "SubAlgo with key " << keys[i] << " not retrieved" ;

    }

    if ( fNuclearFFs.count( ff -> NucleusPDG() ) > 0 ) {

      good_configuration = false ;
      LOG("COHFormFactorMap", pERROR ) << "Attempt to add a second DeVries form factor for PDG " << ff -> NucleusPDG() ;
    }

    fNuclearFFs[ ff -> NucleusPDG() ] = ff ;

  }  // loop over subalgo

  if ( ! good_configuration ) {
    LOG("COHFormFactorMap", pFATAL ) << "Configuration not good, exiting" ;
    exit ( 78 ) ;
  }

}
//____________________________________________________________________________
