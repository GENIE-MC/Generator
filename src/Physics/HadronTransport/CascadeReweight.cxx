//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 
 Julia Tena-Vidal <j.tena-vidal \at liverpool.ac.uk>
 University of Liverpool 

*/
//____________________________________________________________________________

#include <cstdlib>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Physics/HadronTransport/CascadeReweight.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

#include "Physics/HadronTransport/INukeHadroFates.h" 
#include "Framework/Utils/StringUtils.h" 
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//___________________________________________________________________________
CascadeReweight::CascadeReweight() :
EventRecordVisitorI("genie::CascadeReweight")
{

}
//___________________________________________________________________________
CascadeReweight::CascadeReweight(string config) :
EventRecordVisitorI("genie::CascadeReweight", config)
{

}
//___________________________________________________________________________
CascadeReweight::~CascadeReweight()
{

}
//___________________________________________________________________________
void CascadeReweight::ProcessEventRecord(GHepRecord * evrec) const
{
  if( !evrec ){
    LOG("CascadeReweight", pERROR) << "** Null input!" ;
    return ; 
  }
  // Get Associated weight
  double weight = GetEventWeight( * evrec ) ;
  // Set weight 
  evrec->SetWeight( weight ) ;

  return ;
}
//___________________________________________________________________________
double CascadeReweight::GetEventWeight ( const GHepRecord & event ) const{

  GHepParticle * p = 0;
  TIter event_iter(&event);
  double total_weight = fDefaultWeight ;
  while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
    {  
      if( p->Status() != kIStStableFinalState ) continue;
      // Get particle fate
      int fate = event.Particle(p->FirstMother())->RescatterCode() ;
      const auto map_it = fMapFateWeights.find(fate) ; 

      // Get weight given a pdg code.
      if( map_it != fMapFateWeights.end() ) {
	int pdg_target = p->Pdg() ;
	const auto weight_it = (map_it->second).find(pdg_target) ; 
	if( weight_it != (map_it->second).end() ) {
	  total_weight *= weight_it->second ; 
	} else { 
	  total_weight *= fDefaultWeight ;
	}
      } else { 
	// Fate does not exist. Use default weight
	total_weight *= fDefaultWeight ;
      }
    }// end loop over particles
  return total_weight ; 
}
//___________________________________________________________________________
void CascadeReweight::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void CascadeReweight::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//___________________________________________________________________________
void CascadeReweight::LoadConfig(void)
{
  bool good_config = true ; 

  // Get default weight 
  if( GetConfig().Exists("CascadeReweight-Default-Weight") ) {
    GetParam( "CascadeReweight-Default-Weight", fDefaultWeight ) ;
  } else {
    good_config = false ; 
    LOG("CascadeReweight", pERROR) << "Default weight is not specified " ;
  }

  // Clean maps
  fMapFateWeights.clear();
  // Create vector with list of possible keys (follows the order of the fates enumeration)
  std::map<int,string> map_keys { {kIHNFtUndefined,"Undefined"}, {kIHNFtNoInteraction,"NoInteraction"}, 
				  {kIHNFtCEx,"CEx"}, {kIHNFtElas,"Elastic"}, {kIHNFtInelas,"Inelastic"},
				  {kIHNFtAbs,"Abs"}, {kIHNFtCmp,"Cmp"} } ;


  for ( map<int,string>::iterator it_keys = map_keys.begin(); it_keys != map_keys.end(); it_keys++) {
    std::string to_find = "CascadeReweight-Weight-"+(it_keys->second)+"@Pdg=" ;
    auto kpdg_list = GetConfig().FindKeys( to_find.c_str() ) ;
    for( auto kiter = kpdg_list.begin(); kiter != kpdg_list.end(); ++kiter ) {
      const RgKey & key = *kiter ;
      vector<string> kv = genie::utils::str::Split(key,"=");
      assert(kv.size()==2);
      int pdg_target = stoi( kv[1] );
      if( ! PDGLibrary::Instance()->Find(pdg_target) ) {
	LOG("CascadeReweight", pERROR) << "The target Pdg code associated to " << to_find << " is not valid : " << pdg_target ; 
	good_config = false ; 
	continue ; 
      }
      double weight ;
      GetParam( key, weight ) ;
      std::map<int,double> MapWeight ; 
      MapWeight.insert( std::pair<int,double>( pdg_target, weight ) ) ;
      fMapFateWeights[it_keys->first] = MapWeight ; 
    }
  }  

  if( ! good_config ) {
    LOG("CascadeReweight", pERROR) << "Configuration has failed.";
    exit(78) ;
  }
  
}
//___________________________________________________________________________
