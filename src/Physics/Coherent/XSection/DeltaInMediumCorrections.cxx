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
  

  //  LOG("DeltaInMediumCorrections", pINFO) << "Loaded " << fFBCs.size() << " coeffictients for nucleus " << fPDG ; 
  
}
//____________________________________________________________________________
