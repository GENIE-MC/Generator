//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Corey Reed <cjreed \at nikhef.nl>
         Nikhef - February 4, 2010

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include "IBDXSecMap.h"

#include "Algorithm/AlgConfigPool.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

#include <sstream>

using namespace genie;

ClassImp(IBDXSecMap)

//____________________________________________________________________________
IBDXSecMap::IBDXSecMap() :
   XSecAlgorithmI("genie::IBDXSecMap"),
   fDefaultModel(0)
{

}
//____________________________________________________________________________
IBDXSecMap::IBDXSecMap(string config) :
   XSecAlgorithmI("genie::IBDXSecMap", config),
   fDefaultModel(0)
{

}
//____________________________________________________________________________
IBDXSecMap::~IBDXSecMap()
{

}

//____________________________________________________________________________
void IBDXSecMap::Configure(const Registry & config)
{
   Algorithm::Configure(config);
   this->LoadConfig();
}
//____________________________________________________________________________
void IBDXSecMap::Configure(string config)
{
   Algorithm::Configure(config);
   this->LoadConfig();
}
//____________________________________________________________________________
void IBDXSecMap::LoadConfig(void)
{
   // build the default xsec model according to the options contained
   // in IBDXSecMap.xml and/or UserPhysicsOptions.xml
   
   AlgConfigPool * confp = AlgConfigPool::Instance();
   const Registry * gc = confp->GlobalParameterList();

   // load default global model (should work for all nuclei)
   RgAlg dgmodel =
      fConfig->GetAlgDef("IBDNucXsecModel", gc->GetAlg("IBDNucXsecModel"));
   LOG("IBD", pINFO)
      << "Default IBD cross section model: " << dgmodel;
   fDefaultModel =
      dynamic_cast<const XSecAlgorithmI*>(this->SubAlg("IBDNucXsecModel"));
   assert(fDefaultModel!=0);
   
   // check whether to map according to specific isotopes
   fIsotopesUseSameModel = 
      fConfig->GetBoolDef("IsotopesUseSameModel",
			  gc->GetBool("IsotopesUseSameModel"));
   
   // load refined models for specific nuclei
   for(int Z=1; Z<140; Z++) {
      for(int A=Z; A<3*Z; A++) {
	 std::ostringstream key;
	 const int nucpdg = pdg::IonPdgCode(A,Z);
	 key << "IBDNucXsecModel@Pdg=" << nucpdg;
	 RgKey rgkey = key.str();
	 if (this->GetConfig().Exists(rgkey) || gc->Exists(rgkey)) {
	    RgAlg rgmodel = fConfig->GetAlgDef(rgkey, gc->GetAlg(rgkey));
	    LOG("IBD", pNOTICE)
	       << "Nucleus =" << nucpdg
	       << " -> refined nuclear model: " << rgmodel;
	    const XSecAlgorithmI* model =
	       dynamic_cast<const XSecAlgorithmI*>(this->SubAlg(rgkey));
	    assert(model);
	    const int mapkeyval = (fIsotopesUseSameModel) ? Z : nucpdg;
	    fRefinedModels.
	       insert(map<int,const XSecAlgorithmI*>::value_type(mapkeyval,
								 model));
	}
     }
  }
}
//____________________________________________________________________________
const XSecAlgorithmI* IBDXSecMap::SelectModel(const Target & t) const
{
   // search the map for the PDG code of the target
   // if a refined model is found, return it
   // otherwise return the default xsec model
   
   map<int,const XSecAlgorithmI*>::const_iterator it = fRefinedModels.find(t.Pdg());
   
   if(it != fRefinedModels.end()) return it->second;
   else return fDefaultModel;
}
//____________________________________________________________________________
