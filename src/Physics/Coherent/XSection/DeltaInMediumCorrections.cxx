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
#include "Framework/Messenger/Messenger.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"

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
  
  double kf = fKFTable->FindClosestKF(nucleus_pdgc, nucleon_pdgc);
  return kf ;

}
//____________________________________________________________________________
double DeltaInMediumCorrections::AverageDensity( int nucleus_pdg, int nucleon_pdg ) const {
  
  double rho = std::pow( FermiMomentum( int nucleus_pdg, int nucleon_pdg ), 3 ) * 2. / ( 3 * constants::kPi2 );
  return rho ;

}
//____________________________________________________________________________
double DeltaInMediumCorrections::AverageDensity( int nucleus_pdg ) const {
  
  // this is the average density between proton and neutron matters

  std::array<int,2> pdgs = {kPdgProton, kPdgNeutron} ;
  double rho = 0. ;
  for ( auto pdg : pdgs ) {
    rho += AverageDensity( nucleus_pdg, pdg ) ;
  }
  
  rho /= pdgs.size() ;

  return rho ;
  
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

  //  LOG("DeltaInMediumCorrections", pINFO) << "Loaded " << fFBCs.size() << " coeffictients for nucleus " << fPDG ; 
  
}
//____________________________________________________________________________
