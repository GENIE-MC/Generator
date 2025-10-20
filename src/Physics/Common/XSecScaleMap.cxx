//_________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Original code contributed by J.Tena and M.Roda
 Substantial code refactorizations by the core GENIE group.
*/
//_________________________________________________________________________

#include "Physics/Common/XSecScaleMap.h"
#include "Framework/Utils/StringUtils.h" 
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//_________________________________________________________________________
XSecScaleMap::XSecScaleMap() : 
  XSecScaleI("genie::XSecScaleMap") 
{
  
}
//_________________________________________________________________________
XSecScaleMap::XSecScaleMap(string config) : 
  XSecScaleI("genie::XSecScaleMap",config) 
{
  
}
//_________________________________________________________________________
XSecScaleMap::~XSecScaleMap()
{

}
//_________________________________________________________________________
double XSecScaleMap::GetScaling( const Interaction & interaction ) const {
  // This function accesses the requested Algoritm given the Pdg code and 
  // it retrieves the appropiate scaling.
  // Get Target pdg
  int pdg_target = interaction.InitState().Tgt().Pdg() ;

  const auto it = fXSecScaleMap.find(pdg_target) ;
  if ( it != fXSecScaleMap.end() ) {
    return (it -> second)->GetScaling( interaction ) ;
  }
  if ( fXSecScaleDefault ) {
    // return default 
    return fXSecScaleDefault->GetScaling( interaction ) ;
  }
  return 1. ; 

}

//_________________________________________________________________________

void XSecScaleMap::LoadConfig(void)
{
  bool good_config = true ; 
  fXSecScaleDefault = nullptr ; 
  fXSecScaleMap.clear() ; 

  // Store default value
  static RgKey default_algo_name = "XSecScaleDefaultAlg" ;
  if( GetConfig().Exists(default_algo_name) ) {
    fXSecScaleDefault = dynamic_cast<const XSecScaleI *> ( this->SubAlg(default_algo_name) );
    if( !fXSecScaleDefault ) {
      good_config = false ; 
      LOG("XSecScaleMap", pERROR) << "The subalgorithm with ID  couldn't be casted " ; 
      // << SubAlg(default_algo_name)->Id() << " couldn't be casted " ;
    }  

  }
   
  // Get possible entries to pdg - shift map 
  auto kpdg_list = GetConfig().FindKeys("XSecScaleAlg@Pdg=") ;
  
  for( auto kiter = kpdg_list.begin(); kiter != kpdg_list.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"=");
    assert(kv.size()==2);
    int pdg_target = stoi( kv[1] );
    if( ! PDGLibrary::Instance()->Find(pdg_target) ) {
      good_config = false ; 
      LOG("XSecScaleMap", pERROR) << "The target Pdg code associated is not valid : " << pdg_target ; 
      continue ; 
    }
    
    if( ! pdg::IsIon(pdg_target) ) {
      good_config = false ; 
      LOG("XSecScaleMap", pERROR) << "The target Pdg code does not correspond to a Ion : " << pdg_target ; 
      continue ; 
    } 

    const auto algo = fXSecScaleMap[pdg_target] = dynamic_cast<const XSecScaleI*> ( this->SubAlg( key ) ); 
    if( ! algo ) {
      good_config = false ; 
      LOG("XSecScaleMap", pERROR) << "The subalgorithm " << GetConfig().GetAlg(key).name 
				  << " and target pdg " << pdg_target << " do not exist" ;
      continue ; 
    } 

  }
  
  if( ! good_config ) {
    LOG("XSecScaleMap", pERROR) << "Configuration has failed.";
    exit(78) ;
  }

}

//_________________________________________________________________________
