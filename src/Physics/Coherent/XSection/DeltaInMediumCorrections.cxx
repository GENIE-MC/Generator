//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Marco Roda

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <cmath>
#include <complex>
#include <TMath.h>

#include "Physics/Coherent/XSection/DeltaInMediumCorrections.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"


#include "Framework/Messenger/Messenger.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/ParticleData/PDGCodes.h"

#include "Framework/Utils/StringUtils.h" 

using namespace genie;


//____________________________________________________________________________
DeltaInMediumCorrections::DeltaInMediumCorrections() : 
  Algorithm("genie::DeltaInMediumCorrections")
{

}
//____________________________________________________________________________
DeltaInMediumCorrections::DeltaInMediumCorrections(string config) :
  Algorithm("genie::DeltaInMediumCorrections", config)
{

}
//____________________________________________________________________________
DeltaInMediumCorrections::~DeltaInMediumCorrections()
{

}
//____________________________________________________________________________
double DeltaInMediumCorrections::FermiMomentum( int nucleus_pdg, int nucleon_pdg ) const {
  
  double kf = fKFTable->FindClosestKF(nucleus_pdg, nucleon_pdg);
  return kf ;

}
double DeltaInMediumCorrections::FermiMomentum( int nucleus_pdg ) const {
  // Invert avg. nucleus density to obtain avg. Fermi momentum for consistency.
  return pow( ( 3*constants::kPi2*AverageDensity( nucleus_pdg ) / 2 ), 1./3.) ;
}
//____________________________________________________________________________
double DeltaInMediumCorrections::AverageDensity( int nucleus_pdg, int nucleon_pdg ) const {
  
  double rho = std::pow( FermiMomentum( nucleus_pdg, nucleon_pdg ), 3 ) * 2. / ( 3 * constants::kPi2 );
  return rho ;

}
//____________________________________________________________________________
double DeltaInMediumCorrections::AverageDensity( int nucleus_pdg ) const {

  if ( fKFTable ) { 
    // this is the average density between proton and neutron matters
    std::array<int,2> pdgs = { kPdgProton, kPdgNeutron} ;
    double rho = 0. ;
    for ( auto pdg : pdgs ) {
      rho += AverageDensity( nucleus_pdg, pdg ) ;
    }
    
    rho /= pdgs.size() ;

    return rho ;
  }
  else {
    
    static double __rho = 3. / (4. * constants::kPi * pow( 1.2 * units::fermi, 3 ) ) ; 

    return __rho ;
  }

 
}
//____________________________________________________________________________
double DeltaInMediumCorrections::Sigma( int nucleus_pdg ) const {

  double sigma = -0.5 * fDeltaV0 * AverageDensity( nucleus_pdg ) / fRho0 ; 
  return sigma ; 
}
//____________________________________________________________________________
double DeltaInMediumCorrections::Gamma_vacuum( double p2 ) const {

  double Gamma  = 0.0 ;
  
  if ( p2 > pow( constants::kNucleonMass + constants::kPionMass, 2 ) ) { 
      
    double q2cm = Q2_cm( p2, { constants::kNucleonMass2, constants::kPionMass2 } ) ; 
    
    Gamma = ( fDeltaNCoupling2 * constants::kNucleonMass * pow(q2cm, 3./2 ) ) 
      / ( 6.0 * constants::kPi * constants::kPionMass2 * sqrt(p2) ) ; 

  }
    
  return Gamma;
}
//____________________________________________________________________________
double DeltaInMediumCorrections::I_series( double q ) {

  double I = 1.0;

  if (q != 0) {
    if (q > 1.0) I += -2.0 / (5.0 * q * q) + 9.0 / ( 35.0 * pow(q, 4) ) - 2.0 / ( 21.0 * pow(q, 6) ) ;
    else if (q < 1.0) I += 34.0 / 35.0 * q - 22.0 / 105.0 * q * q * q - 1.0 ;
  }

  return I;
}
//____________________________________________________________________________
double DeltaInMediumCorrections::Gamma_tilde( double p2, int nucleus_pdg ) const {

  double q2cm = Q2_cm( p2, { constants::kNucleonMass2, constants::kPionMass2 } ) ;

  if ( q2cm <= 0. ) return 0. ;
  // somtimes the p is smaller than ( m_N + m_pi ) so the q2cm is negative. 
  // In that case the gamma_vacuum will be 0 anyway, so not point in propagating a NaN in the code
  
  double q_tilde = sqrt( q2cm ) / FermiMomentum(nucleus_pdg) ;

  return Gamma_vacuum(p2) * I_series(q_tilde);
}
//____________________________________________________________________________
std::complex<double> DeltaInMediumCorrections::AverageDirectPropagator( double p2, int nucleus_pdg ) const {

  double mDelta = utils::res::Mass( kP33_1232 ) ;
  double mDelta2 = pow( mDelta, 2 ) ;

  // Simplified form since Sigma.real = 0
  return 1.0 / ( p2 - mDelta2 + std::complex<double>(0,1)*mDelta*( Gamma_tilde(p2, nucleus_pdg) - 2*Sigma(nucleus_pdg)) ) ;
}
//____________________________________________________________________________
std::complex<double> DeltaInMediumCorrections::AverageCrossPropagator( double p2 ) const {

  double mDelta = utils::res::Mass( kP33_1232 ) ;
  double mDelta2 = pow( mDelta, 2 ) ;

  return 1.0 / ( p2 - mDelta2 + std::complex<double>(0,1) * mDelta * Gamma_vacuum(p2) ) ;
}
//____________________________________________________________________________
void DeltaInMediumCorrections::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeltaInMediumCorrections::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeltaInMediumCorrections::LoadConfig(void)
{

  fKFTable = nullptr ;

  // get the Fermi momentum table for relativistic Fermi gas
  string table_name ; 
  if ( GetParam( "FermiMomentumTable", table_name, false ) ) {
    FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
    fKFTable = kftp->GetTable( table_name );
    assert(fKFTable);
    LOG( "DeltaInMediumCorrections", pINFO ) << "Configured using Fermi Momentum table called " << table_name ;
  }

  GetParam( "NCG-Delta-V0", fDeltaV0 ) ; 

  GetParam( "NCG-Rho0", fRho0 ) ;
  // Rho0 is in fm-3, we use it in GeV^3 in the code
  fRho0 /= units::fermi3 ;

  double deltaN_coupling = 0. ;
  GetParam( "Delta-N-Coupling", deltaN_coupling ) ;
  fDeltaNCoupling2 = deltaN_coupling * deltaN_coupling ;

}
//____________________________________________________________________________
double DeltaInMediumCorrections::Q2_cm( double s, 
					const std::array<double,2> & masses2 ) {
  
  const double & m1 = masses2[0] ;
  const double & m2 = masses2[1] ;
  return  ( s*s + pow(m1,2) + pow(m2,2) - 2.0*s*m1 - 2.0*m1*m2 - 2.0*s*m2) / (4.0 * s ) ;

}
