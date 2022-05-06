//----------------------------------------------------------------------------
// Implementation file for NHLStepper.h
// Author: John Plows <komninos-john.plows@physics.ox.ac.uk>
// 10/Dec/2021
//----------------------------------------------------------------------------

#include "Physics/NeutralHeavyLepton/NHLStepper.h"

using namespace genie;
using namespace genie::NHL;

std::map< NHLDecayMode_t, double > NHLSelector::fCounterMap;

// empty the vector of counters
void NHLSelector::emptyCounters( ){
    while( fCounterMap.size() > 0 ){ fCounterMap.erase( fCounterMap.begin() ); }
}

// populate the vector of counters
void NHLSelector::fillCounters( std::map< NHLDecayMode_t, double > gammaMap ){
  if( fCounterMap.size() > 0 ) emptyCounters( );
  RandomGen * Rng = RandomGen::Instance();
  
  for( std::map< NHLDecayMode_t, double >::iterator it = gammaMap.begin(); it != gammaMap.end(); ++it )
  { 
    fCounterMap.insert( std::pair< NHLDecayMode_t, double > ( (*it).first, Rng->RndGen().Uniform() ) ); 
  }
}

// inclusive method
std::vector< double > NHLSelector::PropagateTilDecay( SimpleNHL sh ){

    const std::vector< double > V0X = sh.GetOrigin4VX();
    const std::vector< double > V0P = sh.Get4VP();

    const double COMlifetime = sh.GetCoMLifetime( );
    const double LABlifetime = sh.GetLifetime( );

    LOG( "NHL", pDEBUG )
      << "\n COM, LAB lifetimes = " << COMlifetime << ", " << LABlifetime << " [GeV^{-1}]"
      << "\n COM, LAB lifetimes = " << COMlifetime / genie::units::GeV / genie::units::ns 
      << "," << LABlifetime / genie::units::GeV / genie::units::ns << " [ns]";
    
    // Get seed
    RandomGen * Rng = RandomGen::Instance();
    double ranthrow = Rng->RndGen().Uniform( );

    // Exponential decay <-> N = floor( tau * ln( seed^(-1) ) / dt )
    const double t1 = std::floor( LABlifetime *
				  std::log( 1.0 / ranthrow ) /
				  fdt ) * fdt / genie::units::GeV / genie::units::ns;
    // ensure HNL is pointing the right way
    sh.Set4Momentum( V0P );
    const std::vector< double > betaVec  = sh.GetBetaVec( );

    const double nsTocm = ( genie::units::ns * genie::constants::kLightSpeed ) / genie::units::cm;
    std::vector< double > decay4V; decay4V.emplace_back( V0X.at(0) + t1 );
    decay4V.emplace_back( V0X.at(1) + betaVec.at(0) * t1 * nsTocm );
    decay4V.emplace_back( V0X.at(2) + betaVec.at(1) * t1 * nsTocm );
    decay4V.emplace_back( V0X.at(3) + betaVec.at(2) * t1 * nsTocm );

    LOG( "NHL", pDEBUG )
      << "t1 = " << t1 << " [ns], ranthrow = " << ranthrow << ", betaVec = ( "
      << betaVec.at(0) << ", " << betaVec.at(1) << ", " << betaVec.at(2) << " ) ";

    return decay4V;
}

// Exclusive stepper
void NHLSelector::PropagateAndSelectChannel( SimpleNHL sh ){

    __attribute__((unused)) const double COMlifetime = sh.GetCoMLifetime( );
    const double LABlifetime = sh.GetLifetime( );

    // first get the valid and the interesting channels with their widths
    const std::map< NHLDecayMode_t, double > ValidMap = sh.GetValidChannels( );
    const std::map< NHLDecayMode_t, double > InterestingMap = sh.GetInterestingChannels( );

    // TODO: exception for empty map!
    const int vSize = ValidMap.size();
    if( vSize == 0 ){ LOG( "NHL", pERROR ) << 
	"NHLSelector::PropagateAndSelectChannel:: " <<
	"\nNo valid channels!! This NHL cannot decay, which is bogus." <<
	"\nExiting."; exit( 3 ); }
    
    const int iSize = InterestingMap.size( );
    if( iSize == 0 ){ LOG( "NHL", pERROR ) << 
	"NHLSelector::PropagateAndSelectChannel::" <<
	"\nNo interesting channels!! This NHL is invisible." <<
	"\nExiting."; exit( 3 ); }
    
    // ensure vector of counters is empty, then fill
    emptyCounters( );
    fillCounters( ValidMap ); // we will ignore those counters 
    
    // build map of *when* each valid channel would decay, if it were alone
    std::map< NHLDecayMode_t, double > timeMap, tIntMap;

    for( std::map< NHLDecayMode_t, double >::iterator it = fCounterMap.begin();
	 it != fCounterMap.end(); ++it ){
	// Exponential decay <-> N = floor( tau * ln( seed^(-1) ) / dt )
	const double tmpT = std::floor( LABlifetime *
					std::log( 1.0 / (*it).second ) /
					fdt ) * fdt;
	timeMap.insert( std::pair< NHLDecayMode_t, double >( (*it).first, tmpT ) );

	// only insert in tIntMap if we have an interesting channel!
	if( InterestingMap.find( (*it).first ) != InterestingMap.end() ){
	    tIntMap.insert( std::pair< NHLDecayMode_t, double >( (*it).first, tmpT ) ); }
    }

    // first channel wins.
    double t1 = 99999.99, t2 = 99999.99; NHLDecayMode_t c1;
    for( std::map< NHLDecayMode_t, double >::iterator it = timeMap.begin();
	 it != timeMap.end(); ++it ){
	if( (*it).second < t1 ){ t1 = (*it).second; } } // t1 from all valid times
    for( std::map< NHLDecayMode_t, double >::iterator it = tIntMap.begin();
	 it != tIntMap.end(); ++it ){
	if( (*it).second < t2 ){ t2 = (*it).second; c1 = (*it).first; } } // c1 interesting

    const std::vector< double > betaVec  = sh.GetBetaVec( );
    const std::vector< double > origin4V = sh.GetOrigin4VX( );
    std::vector< double > decay4V; decay4V.emplace_back( origin4V.at(0) + t1 );
    decay4V.emplace_back( origin4V.at(1) + betaVec.at(1) * t1 * genie::constants::kLightSpeed );
    decay4V.emplace_back( origin4V.at(2) + betaVec.at(2) * t1 * genie::constants::kLightSpeed );
    decay4V.emplace_back( origin4V.at(3) + betaVec.at(3) * t1 * genie::constants::kLightSpeed );

    sh.SetDecayVtx( decay4V ); sh.SetDecayMode( c1 );
}
