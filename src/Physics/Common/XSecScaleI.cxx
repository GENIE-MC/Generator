//_________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Original code contributed by J.Tena and M.Roda
 Substantial code refactorizations by the core GENIE group.
*/
//_________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Common/XSecScaleI.h"
#include "Framework/Utils/StringUtils.h" 
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//_________________________________________________________________________
XSecScaleI::XSecScaleI() : 
  Algorithm("genie::XSecScaleI") 
{
  
}
//_________________________________________________________________________
XSecScaleI::XSecScaleI(string config) : 
  Algorithm("genie::XSecScaleI",config) 
{
  
}
//_________________________________________________________________________
XSecScaleI::~XSecScaleI()
{

}
//_________________________________________________________________________
void XSecScaleI::Configure(const Registry & config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void XSecScaleI::Configure(string config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//_________________________________________________________________________
double XSecScaleI::GetScaling( const Interaction & interaction ) const {
  // This function accesses the Algoritm given the Pdg code and 
  // retrieves the appropiate scaling.
  // Get Target pdg
  int pdg_target = interaction.InitState().Tgt().Pdg() ;
  
  const auto it = fXSecScaleMap.find(pdg_target) ;
  if ( it != fXSecScaleMap.end() ) {
    return (it -> second)->GetScaling( interaction ) ;
  } else {
    // return default 
    return fXSecScaleDefault ;
  }
  
} 
//_________________________________________________________________________

void XSecScaleI::LoadConfig(void)
{
  bool good_config = true ; 

  // Store default value
  GetParam( "XSecScaleDefault", fXSecScaleDefault ) ;

  // Clear map
  fXSecScaleMap.clear() ; 

  // Get possible entries to pdg - shift map 
  auto kpdg_list = GetConfig().FindKeys("XSecScaleI@Pdg=") ;

  for( auto kiter = kpdg_list.begin(); kiter != kpdg_list.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"=");
    assert(kv.size()==2);
    int pdg_target = stoi( kv[1] );
    if( ! PDGLibrary::Instance()->Find(pdg_target) ) {
      LOG("XSecScaleI", pERROR) << "The target Pdg code associated is not valid : " << pdg_target ; 
      good_config = false ; 
      continue ; 
    }
    
    if( ! pdg::IsIon(pdg_target) ) {
      LOG("XSecScaleI", pERROR) << "The target Pdg code does not correspond to a Ion : " << pdg_target ; 
      good_config = false ; 
      continue ; 
    } 

    fXSecScaleMap[pdg_target] = (XSecScaleI*) ( this->SubAlg( key ) ); 
    if( ! fXSecScaleMap[pdg_target] ) {
      LOG("XSecScaleI", pERROR) << "The subalgorithm with ID " << fXSecScaleMap[pdg_target]->Id() 
				<< " and target pdg " << pdg_target << " does not exist" ;
      good_config = false ; 
      continue ; 
    } 

  }

  if( ! good_config ) {
    LOG("XSecScaleI", pERROR) << "Configuration has failed.";
    exit(78) ;
  }

}

//_________________________________________________________________________
