//____________________________________________________________________________

#include <algorithm>

#include "ReWeight/GReWeightNuXSecParams.h"

using namespace genie;
using namespace genie::rew;

//____________________________________________________________________________
GReWeightNuXSecParams::GReWeightNuXSecParams()
{
  this->Init();
}
//____________________________________________________________________________
GReWeightNuXSecParams::~GReWeightNuXSecParams()
{

}
//____________________________________________________________________________
double GReWeightNuXSecParams::DefValue(GSyst_t syst) const
{
  map<GSyst_t, double>::const_iterator iter = fDefParams.find(syst);
  if(iter != fDefParams.end()) return iter->second;
  else return 0;
}
//____________________________________________________________________________
double GReWeightNuXSecParams::CurValue(GSyst_t syst) const
{
  map<GSyst_t, double>::const_iterator iter = fCurParams.find(syst);
  if(iter != fCurParams.end()) return iter->second;
  else return 0;
}
//____________________________________________________________________________
double GReWeightNuXSecParams::CurTwkDial(GSyst_t syst) const
{
  map<GSyst_t, double>::const_iterator iter = fCurTwkDial.find(syst);
  if(iter != fCurTwkDial.end()) return iter->second;
  else return 0;
}
//____________________________________________________________________________
bool GReWeightNuXSecParams::IsIncluded(GSyst_t syst) const
{
  map<GSyst_t, bool>::const_iterator iter = fIsTweaked.find(syst);
  if(iter != fIsTweaked.end()) return true;
  else return false;
}
//____________________________________________________________________________
bool GReWeightNuXSecParams::IsTweaked(GSyst_t syst) const
{
  map<GSyst_t, bool>::const_iterator iter = fIsTweaked.find(syst);
  if(iter != fIsTweaked.end()) return iter->second;
  else return false;
}
//____________________________________________________________________________
bool GReWeightNuXSecParams::IsTweaked(void) const
{
  map<GSyst_t, bool>::const_iterator iter = fIsTweaked.begin();
  for( ; iter != fIsTweaked.end(); ++iter)
  {
    if(iter->second) return true;
  }
  return false;
}
//____________________________________________________________________________
void GReWeightNuXSecParams::Reset(GSyst_t syst) 
{
  if(this->IsIncluded(syst)) {
    double def_param = this->DefValue(syst);

    this->SetCurValue    (syst, def_param);
    this->SetCurTwkDial  (syst, 0.);
    this->SetTweakedFlag (syst, false);
  }
}
//____________________________________________________________________________
void GReWeightNuXSecParams::SetDefValue(GSyst_t syst, double value)
{
  fDefParams.insert( map<GSyst_t, double>::value_type(syst, value) );     
}
//____________________________________________________________________________
void GReWeightNuXSecParams::SetCurValue(GSyst_t syst, double value)
{
  fCurParams.insert( map<GSyst_t, double>::value_type(syst, value) );     
}
//____________________________________________________________________________
void GReWeightNuXSecParams::SetCurTwkDial(GSyst_t syst, double value)
{
  fCurTwkDial.insert( map<GSyst_t, double>::value_type(syst, value) );     
}
//____________________________________________________________________________
void GReWeightNuXSecParams::SetTweakedFlag(GSyst_t syst, bool value)
{
  fIsTweaked.insert( map<GSyst_t, bool>::value_type(syst, value) );  
}
//____________________________________________________________________________
void GReWeightNuXSecParams::Init(void)
{
  fDefParams.clear ();
  fCurParams.clear ();
  fCurTwkDial.clear();
  fIsTweaked.clear ();
}      
//____________________________________________________________________________
