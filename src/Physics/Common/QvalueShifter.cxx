//_________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Original code contributed by J.Tena and M.Roda
 Substantial code refactorizations by the core GENIE group.
*/
//_________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Common/QvalueShifter.h"
#include "Framework/Utils/StringUtils.h" 
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//_________________________________________________________________________
QvalueShifter::QvalueShifter() : 
  Algorithm("genie::QvalueShifter") 
{
  
}
//_________________________________________________________________________
QvalueShifter::QvalueShifter(string config) : 
  Algorithm("genie::QvalueShifter",config) 
{
  
}
//_________________________________________________________________________
QvalueShifter::~QvalueShifter()
{

}
//_________________________________________________________________________
void QvalueShifter::Configure(const Registry & config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void QvalueShifter::Configure(string config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//_________________________________________________________________________

  double QvalueShifter::Shift( const Target & target ) const {
    // Get Target pdg
    int pdg_target = target.Pdg() ;

    const auto it = fRelShift.find(pdg_target) ;
    if ( it != fRelShift.end() ) {
      return it -> second ;
    } else {
      // return default 
      return fRelShiftDefault ;
    }

  } 

//_________________________________________________________________________

  double QvalueShifter::Shift( const Interaction & interaction ) const {
    // This function allows the flexibility to add shift as a function of 
    // the interaction type.
    // Right now, we just call the default :
    return Shift( interaction.InitState().Tgt() ) ;    
  } 

//_________________________________________________________________________

void QvalueShifter::LoadConfig(void)
{
  bool good_config = true ;
  // Store default value
  GetParam( "QvalueShiftDefault", fRelShiftDefault ) ;

  // Clear map
  fRelShift.clear() ; 

  // Get possible entries to pdg - shift map 
  auto kpdg_list = GetConfig().FindKeys("QvalueShift@Pdg=") ;

  for( auto kiter = kpdg_list.begin(); kiter != kpdg_list.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"=");
    assert(kv.size()==2);
    int pdg_target = stoi( kv[1] );
    if( ! PDGLibrary::Instance()->Find(pdg_target) ) {
      LOG("QvalueShifter", pERROR) << "The target Pdg code associated to the QvalueShift is not valid : " << pdg_target ; 
      good_config = false ; 
      continue ; 
    }

    if( ! pdg::IsIon(pdg_target) ) {
      LOG("QvalueShifter", pERROR) << "The target Pdg code does not correspond to a Ion : " << pdg_target ; 
      good_config = false ; 
      continue ; 
    } 
   GetParam( key, fRelShift[pdg_target] ) ;
  }

  if( ! good_config ) {
    LOG("QvalueShifter", pERROR) << "Configuration has failed.";
    exit(78) ;
  }

}

//_________________________________________________________________________
