//----------------------------------------------------------------------------
// Implementation file for NHLDecaySelector.h
// Author: John Plows <komninos-john.plows@physics.ox.ac.uk>
// 09/Dec/2021
//----------------------------------------------------------------------------
// TODO: bound checking, consistency with scores (final entry == 1 ?)
//----------------------------------------------------------------------------

#include "Physics/NeutralHeavyLepton/NHLDecaySelector.h"

#include "Framework/EventGen/EVGThreadException.h"

using namespace genie::NHL;

// Takes parameter space, outputs all available channels + widths
std::map< NHLDecayMode_t, double > NHLSelector::GetValidChannelWidths( const double M, const double Ue42, const double Umu42, const double Ut42, const bool IsMajorana ){

    std::map< NHLDecayMode_t, double > allChannels;

    // invisible decay is always possible
    double GINV = 0.0;
    if( fDecayGammas[0] < 0.0 ){
      GINV = NHLSelector::DWidth_Invisible( M, Ue42, Umu42, Ut42 );
      fDecayGammas[0] = GINV;
      LOG("NHL", pDEBUG)
	<< " Invisible decay gamma = " << GINV;
    } else GINV = fDecayGammas[0];
    allChannels.insert( std::pair< NHLDecayMode_t, double >( kNHLDcyNuNuNu, GINV ) );

    assert( GINV >= 0.0 );

    // nu-e-e is next lightest
    if( M < 2.0 * genie::constants::kElectronMass ) return allChannels;

    double GNEE = 0.0;
    if( fDecayGammas[1] < 0.0 ){
      GNEE = NHLSelector::DWidth_SameLepton( M, Ue42, Umu42, Ut42, genie::constants::kElectronMass, false );
      fDecayGammas[1] = GNEE;
      LOG("NHL", pDEBUG)
	<< " Nu-e-e gamma = " << GNEE;
    } else GNEE = fDecayGammas[1];
    allChannels.insert( std::pair< NHLDecayMode_t, double >( kNHLDcyNuEE, GNEE ) );

    assert( GNEE >= 0.0 );

    // nu-e-mu is next lightest
    if( M < genie::constants::kElectronMass + genie::constants::kMuonMass ) return allChannels;

    double GNEM = 0.0;
    if( fDecayGammas[2] < 0.0 ){
      GNEM = NHLSelector::DWidth_DiffLepton( M, Ue42, Umu42, IsMajorana ); 
      fDecayGammas[2] = GNEM;
      LOG("NHL", pDEBUG)
	<< " Nu-e-mu gamma = " << GNEM;
    } else GNEM = fDecayGammas[2];
    allChannels.insert( std::pair< NHLDecayMode_t, double >( kNHLDcyNuMuE, GNEM ) );

    assert( GNEM >= 0.0 );

    // pi0-nu is next lightest
    if( M < genie::constants::kPi0Mass ) return allChannels;

    double GP0N = 0.0;
    if( fDecayGammas[3] < 0.0 ){
      GP0N = NHLSelector::DWidth_PiZeroAndNu( M, Ue42, Umu42, Ut42 );
      fDecayGammas[3] = GP0N;
      LOG("NHL", pDEBUG)
	<< " Pi0-nu gamma = " << GP0N;
    } else GP0N = fDecayGammas[3];
    allChannels.insert( std::pair< NHLDecayMode_t, double >( kNHLDcyPi0Nu, GP0N ) );

    assert( GP0N >= 0.0 );

    // pi-e is next lightest
    if( M < genie::constants::kPionMass + genie::constants::kElectronMass ) return allChannels;

    double GPIE = 0.0;
    if( fDecayGammas[4] < 0.0 ){
      GPIE = NHLSelector::DWidth_PiAndLepton( M, Ue42, genie::constants::kElectronMass );
      fDecayGammas[4] = GPIE;
      LOG("NHL", pDEBUG)
	<< " Pi-e gamma = " << GPIE;
    } else GPIE = fDecayGammas[4];
    allChannels.insert( std::pair< NHLDecayMode_t, double >( kNHLDcyPiE, GPIE) );

    assert( GPIE >= 0.0 );

    // nu-mu-mu is next lightest
    if( M < 2.0 * genie::constants::kMuonMass ) return allChannels;

    double GNMM = 0.0;
    if( fDecayGammas[5] < 0.0 ){
      GNMM = NHLSelector::DWidth_SameLepton( M, Ue42, Umu42, Ut42, genie::constants::kMuonMass, false );
      fDecayGammas[5] = GNMM;
      LOG("NHL", pDEBUG)
	<< " Nu-mu-mu gamma = " << GNMM;
    } else GNMM = fDecayGammas[5];
    allChannels.insert( std::pair< NHLDecayMode_t, double >( kNHLDcyNuMuMu, GNMM ) );

    assert( GNMM >= 0.0 );

    // pi-mu is next lightest
    if( M < genie::constants::kPionMass + genie::constants::kMuonMass ) return allChannels;

    double GPIM = 0.0;
    if( fDecayGammas[6] < 0.0 ){
      GPIM = NHLSelector::DWidth_PiAndLepton( M, Umu42, genie::constants::kMuonMass );
      fDecayGammas[6] = GPIM;
      LOG("NHL", pDEBUG)
	<< " Pi-mu gamma  = " << GPIM;
    } else GPIM = fDecayGammas[6];
    allChannels.insert( std::pair< NHLDecayMode_t, double >( kNHLDcyPiMu, GPIM ) );

    assert( GPIM >= 0.0 );

    // pi0-pi0-nu is next lightest
    if( M < 2.0 * genie::constants::kPi0Mass ) return allChannels;

    double GP02 = 0.0;
    if( fDecayGammas[7] < 0.0 ){
      GP02 = NHLSelector::DWidth_Pi0Pi0Nu( M, Ue42, Umu42, Ut42 );
      fDecayGammas[7] = GP02;
      LOG("NHL", pDEBUG)
	<< " Pi0-pi0-nu gamma = " << GP02;
    } else fDecayGammas[7] = GP02;
    allChannels.insert( std::pair< NHLDecayMode_t, double >( kNHLDcyPi0Pi0Nu, GP02 ) );

    assert( GP02 >= 0.0 );

    // pi-pi0-e is next lightest
    if( M < genie::constants::kPionMass + genie::constants::kPi0Mass + genie::constants::kElectronMass ) return allChannels;

    double GP0E = 0.0;
    if( fDecayGammas[8] < 0.0 ){
      GP0E = NHLSelector::DWidth_PiPi0Ell( M, genie::constants::kElectronMass, Ue42, Umu42, Ut42, true );
      fDecayGammas[8] = GP0E;
      LOG("NHL", pDEBUG)
	<< " Pi-pi0-e gamma = " << GP0E;
    } else GP0E = fDecayGammas[8];
    allChannels.insert( std::pair< NHLDecayMode_t, double >( kNHLDcyPiPi0E, GP0E ) );

    assert( GP0E >= 0.0 );

    // pi-pi0-mu is next lightest
    if( M < genie::constants::kPionMass + genie::constants::kPi0Mass + genie::constants::kMuonMass ) return allChannels;

    double GP0M = 0.0;
    if( fDecayGammas[9] < 0.0 ){
      GP0M = NHLSelector::DWidth_PiPi0Ell( M, genie::constants::kMuonMass, Ue42, Umu42, Ut42, false );
      fDecayGammas[9] = GP0M;
      LOG("NHL", pDEBUG)
	<< " Pi-pi0-mu gamma = " << GP0M;
    } else GP0M = fDecayGammas[9];
    allChannels.insert( std::pair< NHLDecayMode_t, double >( kNHLDcyPiPi0Mu, GP0M ) );

    assert( GP0M >= 0.0 );

    //all done! Return
    return allChannels;
}

// Calculates the *total* decay width from all the valid channels
double NHLSelector::GetTotalDecayWidth( std::map< NHLDecayMode_t, double > gammaMap ) {
    
    double totGamma = 0.0;

    for( std::map< NHLDecayMode_t, double >::iterator it = gammaMap.begin(); it != gammaMap.end(); ++it ){ totGamma += (*it).second; }

    LOG("NHL", pDEBUG)
      << " Total gamma from N_channels = " << gammaMap.size()
      << " is = " << totGamma << " [GeV]"
      << " or = " << totGamma * genie::units::GeV * genie::units::ns << " [ns^{-1}]";

    return totGamma;
}

// Returns lifetime of particle with mass and couplings
double NHLSelector::CalcCoMLifetime( const double M, const double Ue42, const double Umu42, const double Ut42, const bool IsMajorana ){

    std::map< NHLDecayMode_t, double > allChannels = NHLSelector::GetValidChannelWidths( M, Ue42, Umu42, Ut42, IsMajorana );
    double totGamma = NHLSelector::GetTotalDecayWidth( allChannels );
    return 1.0 / totGamma; // GeV^{-1}
}

// let's pick the interesting channels
std::map< NHLDecayMode_t, double > NHLSelector::SetInterestingChannels( std::vector< NHLDecayMode_t > intChannels, std::map< NHLDecayMode_t, double > gammaMap ){

    std::map< NHLDecayMode_t, double > interestingMap;

    for( std::vector< NHLDecayMode_t >::iterator it = intChannels.begin(); it != intChannels.end(); ++it ){
	NHLDecayMode_t decType = (*it);
	double gamma = gammaMap.find( (*it) )->second;
	interestingMap.insert( std::pair< NHLDecayMode_t, double >( decType, gamma ) );
    }
    return interestingMap;
} // this is now a reduced map with only the channels we want to decay HNL to

// and transform decay widths to branching ratios (probabilities)
std::map< NHLDecayMode_t, double > NHLSelector::GetProbabilities( std::map< NHLDecayMode_t, double > gammaMap ){

    double totGamma = GetTotalDecayWidth( gammaMap );
    std::map< NHLDecayMode_t, double > Pmap;

    // P = Gamma(channel)/Gamma(tot)
    for( std::map< NHLDecayMode_t, double >::iterator it = gammaMap.begin(); it != gammaMap.end(); ++it ){
	Pmap.insert( std::pair< NHLDecayMode_t, double >( (*it).first, (*it).second / totGamma ) );
    }
    return Pmap;
}

// choose a particular channel to decay to
NHLDecayMode_t NHLSelector::SelectChannelInclusive( std::map< NHLDecayMode_t, double > Pmap, double ranThrow ){

    // in inclusive method, decay is factorised in three parts:
    // a) Decay vertex placement
    // b) Decay product selection: construct 
    // c) Decay product kinematics
    // This method does only b)

    // first get P(all interesting channels)
    double PInt = 0.0, all_before = 0.0;
    NHLDecayMode_t selectedChannel = kNHLDcyNull;
    
    for( std::map< NHLDecayMode_t, double >::iterator it = Pmap.begin(); it != Pmap.end(); ++it ){ PInt += (*it).second; }

    // compare ranThrow to P(channel)/PInt + all_before
    // if all_before + P(channel)/PInt >= ranThrow then select this channel
    // don't break, check if scores add up to 1!
    for( std::map< NHLDecayMode_t, double >::iterator it = Pmap.begin(); it != Pmap.end(); ++it ){
	double modP = (*it).second / PInt;
	if( all_before < ranThrow &&
	    all_before + modP >= ranThrow ) selectedChannel = (*it).first;
	all_before += modP;
    }
    // TODO all_before == 1?

    return selectedChannel;
}
