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
  
  // this is the average density between proton and neutron matters
  std::array<int,2> pdgs = { kPdgProton, kPdgNeutron} ;
  double rho = 0. ;
  for ( auto pdg : pdgs ) {
    rho += AverageDensity( nucleus_pdg, pdg ) ;
  }
  
  rho /= pdgs.size() ;

  return rho ;
  
}
//____________________________________________________________________________
double DeltaInMediumCorrections::Sigma( int nucleus_pdg ) const {

  // TODO verify there isn't an extra factor of 1/2 in Sigma compared to slides
  double sigma = -0.5 * fDeltaV0 * AverageDensity( nucleus_pdg ) / fRho0 ; 
  return sigma ; 
}
//____________________________________________________________________________
double DeltaInMediumCorrections::Gamma_vacuum( double p2 ) const {

  double mn     = constants::kNucleonMass ;
  double mpi    = constants::kPionMass ;

  double p      = sqrt(p2) ;
  double Gamma  = 0.0 ;

  double qcm = sqrt(p2*p2 + pow(mpi,4) + pow(mn,4) - 2.0*p2*mpi*mpi - 2.0*mpi*mpi*mn*mn - 2.0*p2*mn*mn) / (2.0 * p) ;

  if(p2 > (mn + mpi)*(mn + mpi)) {
    Gamma = 1.0 / ( 6.0*constants::kPi ) * ( fDeltaNCoupling/mpi )*( fDeltaNCoupling/mpi )*mn / p*pow(qcm, 3) ;
  }

  return Gamma;
}
//____________________________________________________________________________
double DeltaInMediumCorrections::I_series( double q ) const {

  double I = 1.0;

  if (q != 0) {
    if (q > 1.0) I += -2.0 / (5.0 * q * q) + 9.0 / ( 35.0 * pow(q, 4) ) - 2.0 / ( 21.0 * pow(q, 6) ) ;
    else if (q < 1.0) I += 34.0 / ( 35.0 * q) - 22.0 / (105.0 * q * q * q) - 1.0 ;
  }

  return I;
}
//____________________________________________________________________________
double DeltaInMediumCorrections::Gamma_tilde( double p2, int nucleus_pdg ) const {

  double mn  = constants::kNucleonMass ;
  double mpi = constants::kPionMass ;

  double qcm = sqrt(p2*p2 + pow(mpi,4) + pow(mn,4) - 2.0*p2*mpi*mpi - 2.0*mpi*mpi*mn*mn - 2.0*p2*mn*mn) / ( 2.0 * sqrt(p2) ) ;
  double q_tilde = qcm / FermiMomentum(nucleus_pdg) ;

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

  // get the Fermi momentum table for relativistic Fermi gas
  string table_name ; 
  GetParam( "FermiMomentumTable", table_name ) ;
  
  fKFTable = nullptr ;

  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  fKFTable = kftp->GetTable( table_name );

  assert(fKFTable);

  GetParam( "NCG-Delta-V0", fDeltaV0 ) ; 

  GetParam( "NCG-Rho0", fRho0 ) ;

  GetParam( "Delta-N-Coupling", fDeltaNCoupling ) ;
    
  LOG("DeltaInMediumCorrections", pINFO) << "DeltaV0 " << fDeltaV0 << " Rho0 " << fRho0 << " Delta Coupling " << fDeltaNCoupling;
  
}
//____________________________________________________________________________
