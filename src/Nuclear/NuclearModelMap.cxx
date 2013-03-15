//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

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

  // load refined models for specific nuclei
  for(int Z=1; Z<140; Z++) {
    for(int A=Z; A<3*Z; A++) {
      ostringstream key;
      key << "NuclearModel@Pdg=" << pdg::IonPdgCode(A,Z);
      RgKey rgkey = key.str();
      if (this->GetConfig().Exists(rgkey) || gc->Exists(rgkey)) {
        RgAlg rgmodel = fConfig->GetAlgDef(rgkey, gc->GetAlg(rgkey));
        LOG("Nuclear", pNOTICE) 
          << "Nucleus =" << pdg::IonPdgCode(A,Z) 
                         << " -> refined nuclear model: " << rgmodel;
        const NuclearModelI * model = 
              dynamic_cast<const NuclearModelI *> (this->SubAlg(rgkey));
        assert(model);
        fRefinedModels.insert(map<int,const NuclearModelI*>::value_type(Z,model));
      }
    }
  }
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
