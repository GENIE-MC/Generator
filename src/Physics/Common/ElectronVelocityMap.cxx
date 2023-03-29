///____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 
 \brief  It provides the correct ElectronVelociy model to be used 
         for a given atom. 

  \author   Marco Roda <mroda \at liverpool.ac.uk>
          University of Liverpool
  
  \created March 28, 2023

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include "Physics/Common/ElectronVelocityMap.h"

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

ElectronVelocityMap::ElectronVelocityMap() :
  ElectronVelocity("genie::ElectronVelocityMap", "Default") { ; }
//___________________________________________________________________________
ElectronVelocityMap::ElectronVelocityMap(string config) :
ElectronVelocity("genie::ElectronVelocityMap", config) { ; }
//___________________________________________________________________________
void ElectronVelocityMap::Configure(string config)
{
  Algorithm::Configure(config);
  
  Registry * algos = AlgConfigPool::Instance() -> GlobalParameterList() ;
  Registry r( "ElectronVelocitylMap", false ) ;

  // copy in local pool relevant configurations
  RgIMap entries = algos -> GetItemMap();
  const std::string keyStart = "ElectronVelocity";
  for( RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it ) {

    if( it -> first.compare(0, keyStart.size(), keyStart.c_str()) == 0 ) {
      r.Set( it -> first, algos -> GetAlg(it->first ) ) ;
    }

  }

  this->LoadConfig();
}
//___________________________________________________________________________
void ElectronVelocityMap::LoadConfig() {
  fDefGlobalVelocity = nullptr;
  fSpecificModels.clear();

  // load default global model (should work for all nuclei)
  RgAlg dgmodel ;
  GetParam( "ElectronVelocity", dgmodel ) ;

  LOG("ElcetronVelocity", pINFO)
    << "Default global electron velocity model: " << dgmodel;
  fDefGlobalVelocity = dynamic_cast<const ElectronVelocity *> ( this -> SubAlg( "ElectronVelocity" ) ) ;
  assert(fDefGlobalVelocity);

  // We're looking for keys that match this string
  const std::string keyStart = "ElectronVelocity@Pdg=";
  // Looking in both of these registries
  RgIMap entries = GetConfig().GetItemMap();

  for(RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it){
    const std::string& key = it->first;
    // Does it start with the right string?
    if(key.compare(0, keyStart.size(), keyStart.c_str()) == 0){
      // The rest is the PDG code
      const int pdg = atoi(key.c_str()+keyStart.size());
      const int Z = pdg::IonPdgCodeToZ(pdg);
      //const int A = pdg::IonPdgCodeToA(pdg);

      RgAlg rgmodel = GetConfig().GetAlg(key) ;
      LOG("ElectronVelocity", pNOTICE)
        << "Atom =" << pdg
        << " -> refined velocity model: " << rgmodel;
      const ElectronVelocity * model =
        dynamic_cast<const ElectronVelocity *> (
          this -> SubAlg(key) ) ;
      assert(model);
      fSpecificModels.insert(map<int,const ElectronVelocity*>::value_type(Z,model));
    }
  }
}
//___________________________________________________________________________
void ElectronVelocityMap::InitializeVelocity(Interaction & interaction) const {

  const auto & model = SelectModel(interaction.InitState().Tgt());

  model.InitializeVelocity(interaction);

}
//___________________________________________________________________________
const ElectronVelocity & ElectronVelocityMap::SelectModel(const Target & t) const {
  auto it = fSpecificModels.find(t.Z());

  if ( it != fSpecificModels.end()) return *(it->second) ;
  else return * fDefGlobalVelocity;
}

