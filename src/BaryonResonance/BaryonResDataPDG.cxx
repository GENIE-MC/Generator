//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResDataPDG.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Messenger/Messenger.h"

#include <cstdlib>

using std::endl;

using namespace genie;
using namespace genie::utils::res;

//____________________________________________________________________________
BaryonResDataPDG::BaryonResDataPDG() :
BaryonResDataSetI("genie::BaryonResDataPDG")
{

}
//____________________________________________________________________________
BaryonResDataPDG::BaryonResDataPDG(string config) :
BaryonResDataSetI("genie::BaryonResDataPDG", config)
{

}
//____________________________________________________________________________
BaryonResDataPDG::~BaryonResDataPDG()
{
  fResList.Clear();
  fResIdx.clear();
  fResL.clear();
  fIsD.clear();
  fIsN.clear();
  fResMass.clear();
  fResWidth.clear();
  fResNorm.clear();
}
//____________________________________________________________________________
int BaryonResDataPDG::ResonanceIndex(Resonance_t res) const
{
  map<Resonance_t, int>::const_iterator it = fResIdx.find(res);
  if(it == fResIdx.end()){
    LOG("PDGResData", pFATAL) 
         << "No data for resonance = " << AsString(res);
    abort();
  }
  return it->second;
}
//____________________________________________________________________________
int BaryonResDataPDG::OrbitalAngularMom(Resonance_t res) const
{
  map<Resonance_t, int>::const_iterator it = fResL.find(res);
  if(it == fResL.end()){
    LOG("PDGResData", pFATAL) 
         << "No data for resonance = " << AsString(res);
    abort();
  }
  return it->second;
}
//____________________________________________________________________________
bool BaryonResDataPDG::IsDeltaResonance(Resonance_t res) const
{
  map<Resonance_t, bool>::const_iterator it = fIsD.find(res);
  if(it == fIsD.end()){
    LOG("PDGResData", pFATAL) 
         << "No data for resonance = " << AsString(res);
    abort();
  }
  return it->second;
}
//____________________________________________________________________________
bool BaryonResDataPDG::IsNResonance(Resonance_t res) const
{
  map<Resonance_t, bool>::const_iterator it = fIsN.find(res);
  if(it == fIsN.end()){
    LOG("PDGResData", pFATAL) 
         << "No data for resonance = " << AsString(res);
    abort();
  }
  return it->second;
}
//____________________________________________________________________________
double BaryonResDataPDG::Mass(Resonance_t res) const
{
  map<Resonance_t, double>::const_iterator it = fResMass.find(res);
  if(it == fResMass.end()){
    LOG("PDGResData", pFATAL) 
         << "No data for resonance = " << AsString(res);
    abort();
  }
  return it->second;
}
//____________________________________________________________________________
double BaryonResDataPDG::Width(Resonance_t res) const
{
  map<Resonance_t, double>::const_iterator it = fResWidth.find(res);
  if(it == fResWidth.end()){
    LOG("PDGResData", pFATAL) 
         << "No data for resonance = " << AsString(res);
    abort();
  }
  return it->second;
}
//____________________________________________________________________________
double BaryonResDataPDG::BreitWignerNorm(Resonance_t res) const
{
  map<Resonance_t, double>::const_iterator it = fResNorm.find(res);
  if(it == fResNorm.end()){
    LOG("PDGResData", pFATAL) 
         << "No data for resonance = " << AsString(res);
    abort();
  }
  return it->second;
}
//____________________________________________________________________________
void BaryonResDataPDG::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadResonanceData();
}
//____________________________________________________________________________
void BaryonResDataPDG::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadResonanceData();
}
//____________________________________________________________________________
void BaryonResDataPDG::LoadResonanceData(void)
{
  fResList.Clear();
  fResIdx.clear();
  fResL.clear();
  fIsD.clear();
  fIsN.clear();
  fResMass.clear();
  fResWidth.clear();
  fResNorm.clear();

  assert( fConfig->Exists("ResonanceNameList") );
  string resonanes = fConfig->GetString("ResonanceNameList");

  fResList.DecodeFromNameList(resonanes);

  unsigned int nres = fResList.NResonances();
  for(unsigned int ires = 0; ires < nres; ires++) {

     Resonance_t res = fResList.ResonanceId(ires);
     string res_name = AsString(res);

     string ri_key = res_name + "-ResIndex";
     string rl_key = res_name + "-L";
     string id_key = res_name + "-IsDelta";
     string in_key = res_name + "-IsN";
     string rm_key = res_name + "-Mass";
     string rw_key = res_name + "-Width";
     string bw_key = res_name + "-BreitWignerNorm";

     int    idx   = fConfig -> GetInt    (ri_key);
     int    L     = fConfig -> GetInt    (rl_key);
     bool   isD   = fConfig -> GetBool   (id_key);
     bool   isN   = fConfig -> GetBool   (in_key);
     double mass  = fConfig -> GetDouble (rm_key);
     double width = fConfig -> GetDouble (rw_key);
     double bw    = fConfig -> GetDouble (bw_key);

     fResIdx.   insert (map<Resonance_t, int   >::value_type(res, idx)   );
     fResL.     insert (map<Resonance_t, int   >::value_type(res, L)     );
     fIsD.      insert (map<Resonance_t, bool  >::value_type(res, isD)   );
     fIsN.      insert (map<Resonance_t, bool  >::value_type(res, isN)   );
     fResMass.  insert (map<Resonance_t, double>::value_type(res, mass)  );
     fResWidth. insert (map<Resonance_t, double>::value_type(res, width) );
     fResNorm.  insert (map<Resonance_t, double>::value_type(res, bw)    );
  }
}
//____________________________________________________________________________

