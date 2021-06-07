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
  XSecScaleI(config) 
{
  
}
//_________________________________________________________________________
XSecScaleMap::~XSecScaleMap()
{

}
//_________________________________________________________________________
void XSecScaleMap::Configure(const Registry & config)
{
    XSecScaleI::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void XSecScaleMap::Configure(string config)
{
    XSecScaleI::Configure(config);
    this->LoadConfig();
}
//_________________________________________________________________________
double XSecScaleMap::GetScaling( const Interaction & interaction ) const {
  // This function accesses the Algoritm given the Pdg code and 
  // retrieves the appropiate scaling.
  // Get Target pdg
  int pdg_target = interaction.InitState().Tgt().Pdg() ;
  
  const auto it = fXSecScaleMap.find(pdg_target) ;
  if ( it != fXSecScaleMap.end() ) {
    return (it -> second)->GetScaling( interaction ) ;
  } else if ( fXSecScaleDefault ) {
    // return default 
    return fXSecScaleDefault->GetScaling( interaction ) ;
  }
  return 1. ; 

}

//_________________________________________________________________________

void XSecScaleMap::LoadConfig(void)
{
  bool good_config = true ; 

  // Store default value
  if( GetConfig().Exists("XSecScaleMapAlgDefault") ) {
    fXSecScaleDefault = dynamic_cast<const XSecScaleMap *> ( this->SubAlg("XSecScaleMapAlgDefault") );
    if( !fXSecScaleDefault ) {
      good_config = false ; 
      LOG("XSecScaleMap", pERROR) << "The subalgorithm with ID " << fXSecScaleDefault->Id() << " does not exist " ;
    }  
  } 
  
  // Clear map
  fXSecScaleMap.clear() ; 

  // Get possible entries to pdg - shift map 
  auto kpdg_list = GetConfig().FindKeys("XSecScaleMap@Pdg=") ;

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

    fXSecScaleMap[pdg_target] = (XSecScaleMap*) ( this->SubAlg( key ) ); 
    if( ! fXSecScaleMap[pdg_target] ) {
      good_config = false ; 
      LOG("XSecScaleMap", pERROR) << "The subalgorithm with ID " << fXSecScaleMap[pdg_target]->Id() 
				<< " and target pdg " << pdg_target << " does not exist" ;
      continue ; 
    } 

  }

  if( ! good_config ) {
    LOG("XSecScaleMap", pERROR) << "Configuration has failed.";
    exit(78) ;
  }

}

//_________________________________________________________________________
