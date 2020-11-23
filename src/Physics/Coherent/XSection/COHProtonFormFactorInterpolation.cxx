//____________________________________________________________________________
/*
  Copyright (c) 2003-2020, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Author: Marco Roda
  University of Liverpool
  <mroda \at liverpool.ac.uk>

  For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <vector>
#include <utility>
#include <algorithm>

#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/Coherent/XSection/COHProtonFormFactorInterpolation.h"


using namespace genie;


COHProtonFormFactorInterpolation::COHProtonFormFactorInterpolation() :
  COHFormFactorI("genie::COHProtonFormFactorInterpolation", "Default" ),
  fArchive(), 
  fBaseFF(nullptr),
  fAllowExtrapolation(false)
{

}
//____________________________________________________________________________
COHProtonFormFactorInterpolation::COHProtonFormFactorInterpolation(string config) :
  COHFormFactorI("genie::COHProtonFormFactorInterpolation", config),
  fArchive(), 
  fBaseFF(nullptr),
  fAllowExtrapolation(false)
{

}
//____________________________________________________________________________
COHProtonFormFactorInterpolation::~COHProtonFormFactorInterpolation()
{

}
//____________________________________________________________________________
double COHProtonFormFactorInterpolation::ProtonFF( double Q, int pdg ) const {

  if ( ! HasNucleus(pdg) ) return 0. ;

  if ( fBaseFF -> HasNucleus( pdg ) ) return fBaseFF -> ProtonFF( Q, pdg ) ;

  std::vector<int> neighbours = Neighbours( pdg ) ;

  if ( neighbours.size() == 0 ) return 0. ; 
  // this might happen if the nucleus is too far from any known nucleus 
  // configurations might prevent interpolation or extrapolation
  
  if ( neighbours.size() == 1 ) return fBaseFF ->ProtonFF( Q, neighbours[0] ) ;

  // here we interpolate 
  vector<double> ffs( neighbours.size(), 0. ) ;
  vector<int> zs( neighbours.size(), 0 )  ;

  for ( unsigned int i = 0 ; i < neighbours.size() ; i++ ) {
    ffs[i] = fBaseFF -> ProtonFF( Q, neighbours[i] ) ; 
    zs[i] = pdg::IonPdgCodeToZ( neighbours[i] ) ;
  }

  int z = pdg::IonPdgCodeToZ( pdg ) ;

  return Interpolate( zs, ffs, z ) ; 

}
//____________________________________________________________________________
double COHProtonFormFactorInterpolation::NeutronFF( double Q, int pdg ) const {

  double pff = ProtonFF( Q, pdg ) ; 

  if ( pff == 0. ) return pff ;

  int z = pdg::IonPdgCodeToZ( pdg ) ; 
  
  double scale = ( pdg::IonPdgCodeToA( pdg ) - z ) / (double) z ; 

  return scale * pff ;
  
}
//____________________________________________________________________________
bool COHProtonFormFactorInterpolation::HasNucleus( int pdg ) const {

  if ( fBaseFF -> HasNucleus( pdg ) ) return true ; 
  
  int z =  pdg::IonPdgCodeToZ( pdg ) ;

  if ( fAllowExtrapolation ) {
    if ( z >= fArchive.begin()->first ) return true ;
    else return false;
  }
  else {
    if ( z < fArchive.begin()->first ) return false ;
    else if ( z > fArchive.rbegin()->first ) return false ;
    else return true ;
  }
  
  return true ;
}
//____________________________________________________________________________
genie::Range1D_t COHProtonFormFactorInterpolation::QRange( int pdg ) const { 

  if ( ! HasNucleus(pdg) ) return COHFormFactorI::QRange( pdg ) ;

  if ( fBaseFF -> HasNucleus(pdg) ) return fBaseFF -> QRange( pdg ) ;

  auto neighbours = Neighbours( pdg ) ;

  Range1D_t range = fBaseFF -> QRange( neighbours[0] ) ; 
  // we already tested if the nucleus is ok
  // so the neighbours are not empty and it's always possible to get the first 

  for ( unsigned int i = 1 ; i < neighbours.size() ; ++i ) { 
    auto temp_range = fBaseFF -> QRange( neighbours[i] ) ;
    if ( temp_range.min < range.min )  range.min = temp_range.min ; 
    if ( temp_range.max > range.max )  range.max = temp_range.max ; 
  } 

  return range ;
		      
}
//____________________________________________________________________________
std::vector<int> COHProtonFormFactorInterpolation::Neighbours( int pdg ) const {

  int z = pdg::IonPdgCodeToZ( pdg ) ;
  int n = pdg::IonPdgCodeToA( pdg ) - z ;

  std::vector<int> result ; 

  // first we check if a nucleus with the same z is available 
  // in that case that's what we return
  if ( fArchive.count( z ) > 0 ) {
    // find the closest N
    result.push_back( ClosestIsotope( fArchive.at(z), n ) ) ;
    return result ;
  }

  // here we are in the case where we need to extract 2 nuclei from the archive

  auto res = fArchive.rbegin() ;
  
  // then if z > max_z
  // in that case we take the latest nuclei we at our disposal for extrapolation 
  if ( res -> first < z ) { 
    
    result.push_back( ClosestIsotope( res->second, n ) ) ;
    result.push_back( ClosestIsotope( (++res)->second, n ) ) ;
    
    return result ; 
  }

  
  // otherwise we loop down 

  do {
   res++ ;
  } while ( res -> first > z && res != fArchive.rend() ) ;
 
  // since we are iterating backward
  // now res is point to the closes z we have that is smaller than the z we need

  auto second = res ;
  second-- ;
  
  result.push_back( ClosestIsotope( second->second, n ) ) ;
  result.push_back( ClosestIsotope( res->second, n ) ) ;
  
  return result ; 
}
//____________________________________________________________________________
int COHProtonFormFactorInterpolation::ClosestIsotope( const neutron_map & map , 
						      int n ) const {

  if ( n > map.rbegin()->first ) return map.rbegin()->second ;
  
  if ( n < map.begin()->first ) return map.begin()->second ;

  auto min_it = map.begin() ;
  while ( min_it->first < n ) {
    min_it ++ ;
  }
  
  // at this point the iterator min_it points to the closet after n
  // 
  auto max_it = min_it -- ; 
  
  // we now check which one of the two is che closes

  if ( max_it -> first - n < n - min_it -> first ) {
    // the closest is max 
    return max_it -> second ;
  }
  else {
    return min_it -> second ;
  }
}
//____________________________________________________________________________
double COHProtonFormFactorInterpolation::Interpolate( const vector<int> & zs, 
						      const vector<double> & ffs,
						      int final_z ) const  {

  // linear interpolation
  double ret = ffs[0] + (ffs[1] - ffs[0])*(final_z - zs[0])/(double)(zs[1] - zs[0]) ; 

  return ret ;
  
}

//____________________________________________________________________________
void COHProtonFormFactorInterpolation::LoadConfig(void)
{

  fArchive.clear() ;

  bool good_configuration = true ;

  fBaseFF = dynamic_cast<const genie::COHFormFactorI *>( SubAlg( "BaseCOHFormFactor" ) ) ;

  if ( ! fBaseFF ) {
    good_configuration = false ;
    LOG("COHProtonFormFactorInterpolation", pERROR ) << "Invalid Base Form Factor subalgo" ; 
  }

  // load the archive for fast retrieving during operation
  for ( unsigned int z = 1 ; z < 100 ; ++z ) {  
    
    neutron_map temp ;

    for ( unsigned int n = z/2 ; n <= 2*z ; ++n ) {
      
      int pdg = pdg::IonPdgCode( n+z, z ) ;

      TParticlePDG * pdg_particle = PDGLibrary::Instance() -> Find( pdg ) ; 
      if ( ! pdg_particle ) continue ;
      
      LOG("COHProtonFormFactorInterpolation", pINFO ) << "Adding nucleus " << pdg ; 

      temp[n] = pdg ;

    }

    fArchive[z] = std::move(temp) ;
  }
 
  
  if ( fArchive.size() < 2 ) {
    LOG("COHProtonFormFactorInterpolation", pERROR ) << "Not enough FormFactors inside the Map" ;
    good_configuration = false ;
  }

  GetParamDef( "AllowExtrapolation", fAllowExtrapolation, false ) ;

  if ( ! good_configuration ) {
    LOG("COHProtonFormFactorInterpolation", pFATAL ) << "Configuration not good, exiting" ;
    exit ( 78 ) ;
  }

}
//____________________________________________________________________________
