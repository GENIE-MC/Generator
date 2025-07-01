 //____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Mar 18, 2016- Joe Johnston (SD)
   Update GenerateNucleon() and Prob() to accept a radius as the argument,
   and call the corresponding methods in the nuclear model with a radius.

*/
//____________________________________________________________________________


#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________

bool NuclearModelI::GenerateNucleon(const Target & tgt,
                                    double /*hitNucleonRadius*/) const
  {
    return GenerateNucleon(tgt);
  }

double NuclearModelI::Prob(double p, double w, const Target & tgt,
                           double /*hitNucleonRadius*/) const
  {
    return Prob(p,w,tgt);
  }

//____________________________________________________________________________

double NuclearModelI::FermiMomentum( const Target & t, int nucleon_pdg ) const {

  if ( ! fKFTable ) return 0. ; 

  return fKFTable->FindClosestKF( t.Pdg(), nucleon_pdg);

}

//____________________________________________________________________________

double NuclearModelI::LocalFermiMomentum( const Target & t, int nucleon_pdg, 
					  double /*radius*/ ) const {
  return FermiMomentum( t, nucleon_pdg ) ;

}

//____________________________________________________________________________

void NuclearModelI::LoadConfig() {

  string fermi_table_key = "FermiMomentumTable" ;

  // first try to get the Fermi Momentum table from the specific model configurtaion
  if ( ! GetParam( fermi_table_key, fKFTableName, false ) ) {

    // if that fails, the information should come from the Global Config 
    Registry * algos = AlgConfigPool::Instance() -> GlobalParameterList() ;
    fKFTableName = algos -> GetString( fermi_table_key ) ;
    
  }

  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  fKFTable = kftp->GetTable(fKFTableName);

  // Note that model specifications can abvoid the usage of the table
  // but if this configuration is called it's necessary that the table is set.
  assert(fKFTable);
  
}
