//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include <TSystem.h>
#include <TNtupleD.h>
#include <TGraph2D.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearModelMap.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Numerical/RandomGen.h"

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
bool NuclearModelMap::GenerateNucleon(const Target & target) const
{
  const NuclearModelI * nm = this->SelectModel(target);
  if(!nm) return false;

  bool ok = nm->GenerateNucleon(target);

  fCurrRemovalEnergy = nm->RemovalEnergy();
  TVector3 p = nm->Momentum3();
  fCurrMomentum.SetXYZ(p.Px(), p.Py(), p.Pz());
  fFermiMoverInteractionType = nm->GetFermiMoverInteractionType();
  
  return ok;
}
//____________________________________________________________________________
double NuclearModelMap::Prob(double p, double w, const Target & target) const
{
  const NuclearModelI * nm = this->SelectModel(target);
  if(!nm) return 0;

  return nm->Prob(p,w,target);
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
  this->LoadConfig();
}
//____________________________________________________________________________
void NuclearModelMap::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fDefGlobModel = 0;

  // load default global model (should work for all nuclei)
  //
  RgAlg dgmodel =
    fConfig->GetAlgDef("NuclearModel", gc->GetAlg("NuclearModel"));
  LOG("Nuclear", pINFO)
    << "Default global nuclear model: " << dgmodel;
  fDefGlobModel =
    dynamic_cast<const NuclearModelI *> (this->SubAlg("NuclearModel"));
  assert(fDefGlobModel);

  // We're looking for keys that match this string
  const std::string keyStart = "NuclearModel@Pdg=";
  // Looking in both of these registries
  RgIMap entries = fConfig->GetItemMap();
  RgIMap gcEntries = gc->GetItemMap();
  entries.insert(gcEntries.begin(), gcEntries.end());

  for(RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it){
    const std::string& key = it->first;
    // Does it start with the right string?
    if(key.compare(0, keyStart.size(), keyStart.c_str()) == 0){
      // The rest is the PDG code
      const int pdg = atoi(key.c_str()+keyStart.size());
      const int Z = pdg::IonPdgCodeToZ(pdg);
      //const int A = pdg::IonPdgCodeToA(pdg);

      RgAlg rgmodel = fConfig->GetAlgDef(key, gc->GetAlg(key));
      LOG("Nuclear", pNOTICE)
        << "Nucleus =" << pdg
        << " -> refined nuclear model: " << rgmodel;
      const NuclearModelI * model =
        dynamic_cast<const NuclearModelI *> (this->SubAlg(key));
      assert(model);
      fRefinedModels.insert(map<int,const NuclearModelI*>::value_type(Z,model));
    }
  }
 
  LOG("Nuclear", pDEBUG)
    << "Finished LoadConfig";
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
