 //____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Mar 18, 2016- Joe Johnston (SD)
   Update GenerateNucleon() and Prob() to accept a radius as the argument,
   and call the corresponding methods in the nuclear model with a radius.

*/
//____________________________________________________________________________

#include <sstream>

#include <TSystem.h>
#include <TNtupleD.h>
#include <TGraph2D.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/NuclearModelMap.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/RandomGen.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
NuclearModelMap::NuclearModelMap() :
NuclearModelI("genie::NuclearModelMap")
{

}
//____________________________________________________________________________
NuclearModelMap::NuclearModelMap(string config) :
NuclearModelI("genie::NuclearModelMap", config)
{

}
//____________________________________________________________________________
NuclearModelMap::~NuclearModelMap()
{

}
//____________________________________________________________________________
bool NuclearModelMap::GenerateNucleon(const Target & target,
                                      double hitNucleonRadius) const
{
  const NuclearModelI * nm = this->SelectModel(target);
  if(!nm) return false;

  bool ok = nm->GenerateNucleon(target,hitNucleonRadius);

  fCurrRemovalEnergy = nm->RemovalEnergy();
  const TVector3& p  = nm->Momentum3();
  fCurrMomentum.SetXYZ(p.Px(), p.Py(), p.Pz());
  fFermiMoverInteractionType = nm->GetFermiMoverInteractionType();

  return ok;
}
//____________________________________________________________________________
double NuclearModelMap::Prob(double p, double w, const Target & target,
                             double hitNucRadius) const
{
  const NuclearModelI * nm = this->SelectModel(target);
  if(!nm) return 0;

  return nm->Prob(p,w,target,hitNucRadius);
}
//____________________________________________________________________________
NuclearModel_t NuclearModelMap::ModelType(const Target & target) const
{
  const NuclearModelI * nm = this->SelectModel(target);
  if(!nm) return kNucmUndefined;

  return nm->ModelType(target);
}
//____________________________________________________________________________
void NuclearModelMap::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuclearModelMap::Configure(string config)
{
  Algorithm::Configure(config);

  Registry * algos = AlgConfigPool::Instance() -> GlobalParameterList() ;
  Registry r( "NuclearModelMap", false ) ;

  // copy in local pool relevant configurations
  RgIMap entries = algos -> GetItemMap();
  const std::string keyStart = "NuclearModel";
  for( RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it ) {

    if( it -> first.compare(0, keyStart.size(), keyStart.c_str()) == 0 ) {
      r.Set( it -> first, algos -> GetAlg(it->first ) ) ;
    }

  }

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//____________________________________________________________________________
void NuclearModelMap::LoadConfig(void)
{

  fDefGlobModel = 0;
  // load default global model (should work for all nuclei)
  RgAlg dgmodel ;
  GetParam( "NuclearModel", dgmodel ) ;

  LOG("Nuclear", pINFO)
    << "Default global nuclear model: " << dgmodel;
  fDefGlobModel = dynamic_cast<const NuclearModelI *> ( this -> SubAlg( "NuclearModel" ) ) ;
  assert(fDefGlobModel);

  // We're looking for keys that match this string
  const std::string keyStart = "NuclearModel@Pdg=";
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
      LOG("Nuclear", pNOTICE)
        << "Nucleus =" << pdg
        << " -> refined nuclear model: " << rgmodel;
      const NuclearModelI * model =
        dynamic_cast<const NuclearModelI *> (
          this -> SubAlg(key) ) ;
      assert(model);
      fRefinedModels.insert(map<int,const NuclearModelI*>::value_type(Z,model));
    }
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  for (map<int,const NuclearModelI*>::iterator it = fRefinedModels.begin();
      it != fRefinedModels.end(); ++it) {
    LOG("Nuclear", pDEBUG)
      << "Z = " << (*it).first << "; model = " << (*it).second;
  }
#endif
}

//____________________________________________________________________________
const NuclearModelI * NuclearModelMap::SelectModel(const Target & t) const
{
  int Z = t.Z();

  map<int,const NuclearModelI*>::const_iterator it = fRefinedModels.find(Z);

  if(it != fRefinedModels.end()) return it->second;
  else return fDefGlobModel;
}
//____________________________________________________________________________
