//____________________________________________________________________________
/*
  Copyright (c) 2003-2022, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  
  Author: John Plows <komninos-john.plows \at physics.ox.ac.uk>
          University of Oxford
*/
//____________________________________________________________________________

#include "Physics/BeamHNL/HNLDecaySelector.h"

#include "Framework/EventGen/EVGThreadException.h"

using namespace genie::HNL;

// Takes parameter space, outputs all available channels + widths
std::map< HNLDecayMode_t, double > HNLSelector::GetValidChannelWidths( const double M, const double Ue42, const double Umu42, const double Ut42, const bool IsMajorana ){

  // construct an HNLBRFunctions * object to handle the scalings.
  const Algorithm * algBRFunc = AlgFactory::Instance()->GetAlgorithm("genie::HNL::HNLBRFunctions", "Default");
  const HNLBRFunctions * BRFunc = dynamic_cast< const HNLBRFunctions * >( algBRFunc );
  
  std::map< HNLDecayMode_t, double > allChannels;
  
  // invisible decay is always possible
  double GINV = 0.0;
  if( fDecayGammas[0] < 0.0 ){
    GINV = BRFunc->DWidth_Invisible( M, Ue42, Umu42, Ut42 );
    if( IsMajorana ) GINV *= 2.0;
    fDecayGammas[0] = GINV;
    LOG("HNL", pDEBUG)
      << " Invisible decay gamma = " << fDecayGammas[0];
  } else GINV = fDecayGammas[0];
  allChannels.insert( allChannels.begin(), std::pair< HNLDecayMode_t, double >( kHNLDcyNuNuNu, GINV ) );

  assert( GINV >= 0.0 );

  // nu-e-e is next lightest
  if( M < 2.0 * genie::constants::kElectronMass ) return allChannels;

  double GNEE = 0.0;
  if( fDecayGammas[1] < 0.0 ){
    GNEE = BRFunc->DWidth_SameLepton( M, Ue42, Umu42, Ut42, genie::constants::kElectronMass, false );
    if( IsMajorana ) GNEE *= 2.0;
    fDecayGammas[1] = GNEE;
    LOG("HNL", pDEBUG)
      << " Nu-e-e gamma = " << fDecayGammas[1];
  } else GNEE = fDecayGammas[1];
  allChannels.insert( allChannels.begin(), std::pair< HNLDecayMode_t, double >( kHNLDcyNuEE, GNEE ) );

  assert( GNEE >= 0.0 );

  // nu-e-mu is next lightest
  if( M < genie::constants::kElectronMass + genie::constants::kMuonMass ) return allChannels;

  double GNEM = 0.0;
  if( fDecayGammas[2] < 0.0 ){
    GNEM = BRFunc->DWidth_DiffLepton( M, Ue42, Umu42, IsMajorana ); 
    if( IsMajorana ) GNEM *= 2.0;
    fDecayGammas[2] = GNEM;
    LOG("HNL", pDEBUG)
      << " Nu-e-mu gamma = " << fDecayGammas[2];
  } else GNEM = fDecayGammas[2];
  allChannels.insert( allChannels.begin(), std::pair< HNLDecayMode_t, double >( kHNLDcyNuMuE, GNEM ) );

  assert( GNEM >= 0.0 );

  // pi0-nu is next lightest
  if( M < genie::constants::kPi0Mass ) return allChannels;

  double GP0N = 0.0;
  if( fDecayGammas[3] < 0.0 ){
    GP0N = BRFunc->DWidth_PiZeroAndNu( M, Ue42, Umu42, Ut42 );
    if( IsMajorana ) GP0N *= 2.0;
    fDecayGammas[3] = GP0N;
    LOG("HNL", pDEBUG)
      << " Pi0-nu gamma = " << fDecayGammas[3];
  } else GP0N = fDecayGammas[3];
  allChannels.insert( allChannels.begin(), std::pair< HNLDecayMode_t, double >( kHNLDcyPi0Nu, GP0N ) );

  assert( GP0N >= 0.0 );

  // pi-e is next lightest
  if( M < genie::constants::kPionMass + genie::constants::kElectronMass ) return allChannels;

  double GPIE = 0.0;
  if( fDecayGammas[4] < 0.0 ){
    GPIE = BRFunc->DWidth_PiAndLepton( M, Ue42, genie::constants::kElectronMass );
    if( IsMajorana ) GPIE *= 2.0;
    fDecayGammas[4] = GPIE;
    LOG("HNL", pDEBUG)
      << " Pi-e gamma = " << fDecayGammas[4];
  } else GPIE = fDecayGammas[4];
  allChannels.insert( allChannels.begin(), std::pair< HNLDecayMode_t, double >( kHNLDcyPiE, GPIE) );

  assert( GPIE >= 0.0 );

  // nu-mu-mu is next lightest
  if( M < 2.0 * genie::constants::kMuonMass ) return allChannels;

  double GNMM = 0.0;
  if( fDecayGammas[5] < 0.0 ){
    GNMM = BRFunc->DWidth_SameLepton( M, Ue42, Umu42, Ut42, genie::constants::kMuonMass, false );
    if( IsMajorana ) GNMM *= 2.0;
    fDecayGammas[5] = GNMM;
    LOG("HNL", pDEBUG)
      << " Nu-mu-mu gamma = " << fDecayGammas[5];
  } else GNMM = fDecayGammas[5];
  allChannels.insert( allChannels.begin(), std::pair< HNLDecayMode_t, double >( kHNLDcyNuMuMu, GNMM ) );

  assert( GNMM >= 0.0 );

  // pi-mu is next lightest
  if( M < genie::constants::kPionMass + genie::constants::kMuonMass ) return allChannels; 

  double GPIM = 0.0;
  if( fDecayGammas[6] < 0.0 ){
    GPIM = BRFunc->DWidth_PiAndLepton( M, Umu42, genie::constants::kMuonMass );
    if( IsMajorana ) GPIM *= 2.0;
    fDecayGammas[6] = GPIM;
    LOG("HNL", pDEBUG)
      << " Pi-mu gamma  = " << fDecayGammas[6];
  } else GPIM = fDecayGammas[6];
  allChannels.insert( allChannels.begin(), std::pair< HNLDecayMode_t, double >( kHNLDcyPiMu, GPIM ) );

  assert( GPIM >= 0.0 );

  // pi0-pi0-nu is next lightest
  if( M < 2.0 * genie::constants::kPi0Mass ) return allChannels;

  double GP02 = 0.0;
  if( fDecayGammas[7] < 0.0 ){
    GP02 = BRFunc->DWidth_Pi0Pi0Nu( M, Ue42, Umu42, Ut42 );
    if( IsMajorana ) GP02 *= 2.0;
    fDecayGammas[7] = GP02;
    LOG("HNL", pDEBUG)
      << " Pi0-pi0-nu gamma = " << fDecayGammas[7];
  } else fDecayGammas[7] = GP02;
  allChannels.insert( allChannels.begin(), std::pair< HNLDecayMode_t, double >( kHNLDcyPi0Pi0Nu, GP02 ) );

  assert( GP02 >= 0.0 );

  // pi-pi0-e is next lightest
  if( M < genie::constants::kPionMass + genie::constants::kPi0Mass + genie::constants::kElectronMass ) return allChannels;

  double GP0E = 0.0;
  if( fDecayGammas[8] < 0.0 ){
    GP0E = BRFunc->DWidth_PiPi0Ell( M, genie::constants::kElectronMass, Ue42, Umu42, Ut42, true );
    if( IsMajorana ) GP0E *= 2.0;
    fDecayGammas[8] = GP0E;
    LOG("HNL", pDEBUG)
      << " Pi-pi0-e gamma = " << fDecayGammas[8];
  } else GP0E = fDecayGammas[8];
  allChannels.insert( allChannels.begin(), std::pair< HNLDecayMode_t, double >( kHNLDcyPiPi0E, GP0E ) );

  assert( GP0E >= 0.0 );

  // pi-pi0-mu is next lightest
  if( M < genie::constants::kPionMass + genie::constants::kPi0Mass + genie::constants::kMuonMass ) return allChannels;

  double GP0M = 0.0;
  if( fDecayGammas[9] < 0.0 ){
    GP0M = BRFunc->DWidth_PiPi0Ell( M, genie::constants::kMuonMass, Ue42, Umu42, Ut42, false );
    if( IsMajorana ) GP0M *= 2.0;
    fDecayGammas[9] = GP0M;
    LOG("HNL", pDEBUG)
      << " Pi-pi0-mu gamma = " << fDecayGammas[9];
  } else GP0M = fDecayGammas[9];
  allChannels.insert( allChannels.begin(), std::pair< HNLDecayMode_t, double >( kHNLDcyPiPi0Mu, GP0M ) );

  assert( GP0M >= 0.0 );

  //all done! Return
  return allChannels;
}

// Calculates the *total* decay width from all the valid channels
double HNLSelector::GetTotalDecayWidth( std::map< HNLDecayMode_t, double > gammaMap ) {
    
  double totGamma = 0.0;

  for( std::map< HNLDecayMode_t, double >::iterator it = gammaMap.begin(); it != gammaMap.end(); ++it ){ totGamma += (*it).second; }

  LOG("HNL", pDEBUG)
    << " Total gamma from N_channels = " << gammaMap.size()
    << " is = " << totGamma << " [GeV]"
    << " or = " << totGamma * genie::units::GeV * genie::units::ns << " [ns^{-1}]";

  return totGamma;
}

// Returns lifetime of particle with mass and couplings
double HNLSelector::CalcCoMLifetime( const double M, const double Ue42, const double Umu42, const double Ut42, const bool IsMajorana ){

  std::map< HNLDecayMode_t, double > allChannels = HNLSelector::GetValidChannelWidths( M, Ue42, Umu42, Ut42, IsMajorana );
  double totGamma = HNLSelector::GetTotalDecayWidth( allChannels );
  return 1.0 / totGamma; // GeV^{-1}
}

// let's pick the interesting channels
std::map< HNLDecayMode_t, double > HNLSelector::SetInterestingChannels( std::vector< HNLDecayMode_t > intChannels, std::map< HNLDecayMode_t, double > gammaMap ){

  std::map< HNLDecayMode_t, double > interestingMap;

  for( std::vector< HNLDecayMode_t >::iterator it = intChannels.begin(); it != intChannels.end(); ++it ){
    HNLDecayMode_t decType = (*it);
    double gamma = gammaMap.find( (*it) )->second;
    interestingMap.insert( std::pair< HNLDecayMode_t, double >( decType, gamma ) );
  }
  return interestingMap;
} // this is now a reduced map with only the channels we want to decay HNL to

// and transform decay widths to branching ratios (probabilities)
std::map< HNLDecayMode_t, double > HNLSelector::GetProbabilities( std::map< HNLDecayMode_t, double > gammaMap ){

  double totGamma = GetTotalDecayWidth( gammaMap );
  std::map< HNLDecayMode_t, double > Pmap;

  // P = Gamma(channel)/Gamma(tot)
  for( std::map< HNLDecayMode_t, double >::iterator it = gammaMap.begin(); it != gammaMap.end(); ++it ){
    Pmap.insert( std::pair< HNLDecayMode_t, double >( (*it).first, (*it).second / totGamma ) );
  }
  return Pmap;
}

// choose a particular channel to decay to
HNLDecayMode_t HNLSelector::SelectChannelInclusive( std::map< HNLDecayMode_t, double > Pmap, double ranThrow ){

  // in inclusive method, decay is factorised in three parts:
  // a) Decay vertex placement
  // b) Decay product selection: construct 
  // c) Decay product kinematics
  // This method does only b)

  // first get P(all interesting channels)
  double PInt = 0.0, all_before = 0.0;
  HNLDecayMode_t selectedChannel = kHNLDcyNull;
    
  for( std::map< HNLDecayMode_t, double >::iterator it = Pmap.begin(); it != Pmap.end(); ++it ){ PInt += (*it).second; }

  // compare ranThrow to P(channel)/PInt + all_before
  // if all_before + P(channel)/PInt >= ranThrow then select this channel
  // don't break, check if scores add up to 1!
  for( std::map< HNLDecayMode_t, double >::iterator it = Pmap.begin(); it != Pmap.end(); ++it ){
    double modP = (*it).second / PInt;
    if( all_before < ranThrow &&
	all_before + modP >= ranThrow ) selectedChannel = (*it).first;
    all_before += modP;
  }

  return selectedChannel;
}
