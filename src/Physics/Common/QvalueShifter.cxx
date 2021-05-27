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
#include "Physics/Common/QvalueShifter.h"
#include "Framework/Utils/StringUtils.h" 

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

    // Find Pdg in map:
    if ( fRelShift.find(pdg_target) == fRelShift.end() ) {
      return fRelShift.at(pdg_target);
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

    return Shift( * interaction.InitState().TgtPtr() ) ;
    
  } 
//_________________________________________________________________________

void QvalueShifter::LoadConfig(void)
{

  // Store default value
  GetParam( "QvalueShiftDefault", fRelShiftDefault, 0. ) ;

  // Get possible entries to pdg - shift map 
  RgKeyList kpdg_list = GetConfig().FindKeys("QvalueShift@Pdg=") ;
  RgKeyList::const_iterator kiter = kpdg_list.begin();
  for( ; kiter != kpdg_list.end(); ++kiter ) {
    RgKey key = *kiter ; 
    double rshift = GetConfig().GetDouble(key);
    vector<string> kv = genie::utils::str::Split(key,"=");
    assert(kv.size()==2);
    int pdg_target = atoi(kv[1].c_str());
    fRelShift.insert( std::pair<int,double>( pdg_target , rshift ) );
  }

}

//_________________________________________________________________________
