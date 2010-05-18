//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "ReWeight/GSystSet.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GSystSet::GSystSet()
{

}
//_______________________________________________________________________________________
GSystSet::GSystSet(const GSystSet & syst)
{
  this->Copy(syst);
}
//_______________________________________________________________________________________
GSystSet::~GSystSet()
{
  fSystematics.clear();
}
//_______________________________________________________________________________________
void GSystSet::Include(GSyst_t syst)
{
  if(syst == kNullSystematic) return;

  GSystInfo * syst_info = new GSystInfo;
  fSystematics.insert( map<GSyst_t, GSystInfo*>::value_type(syst, syst_info) );     
}
//_______________________________________________________________________________________
int GSystSet::NIncluded(void) const
{
  return fSystematics.size();
}
//_______________________________________________________________________________________
bool GSystSet::IsIncluded(GSyst_t syst) const
{
  return (fSystematics.find(syst) != fSystematics.end());
}
//_______________________________________________________________________________________
vector<genie::rew::GSyst_t> GSystSet::AllIncluded(void)
{
  vector<GSyst_t> svec;

  map<GSyst_t, GSystInfo*>::const_iterator it = fSystematics.begin();
  for( ; it != fSystematics.end(); ++it) {
    GSyst_t syst = it->first;
    svec.push_back(syst);
  }
  return svec;
}
//_______________________________________________________________________________________
void GSystSet::Remove(GSyst_t syst)
{
  fSystematics.erase(syst);
}
//_______________________________________________________________________________________
double GSystSet::CurValue(GSyst_t syst) const
{
  if ( this->IsIncluded(syst) ) return fSystematics.find(syst)->second->CurValue;
  else return 0.;
}
//_______________________________________________________________________________________
double GSystSet::InitValue(GSyst_t syst) const
{
  if ( this->IsIncluded(syst) ) return fSystematics.find(syst)->second->InitValue;
  else return 0.;
}
//_______________________________________________________________________________________
double GSystSet::MinValue(GSyst_t syst) const
{
  if ( this->IsIncluded(syst) ) return fSystematics.find(syst)->second->MinValue;
  else return 0.;
}
//_______________________________________________________________________________________
double GSystSet::MaxValue(GSyst_t syst) const
{
  if ( this->IsIncluded(syst) ) return fSystematics.find(syst)->second->MaxValue;
  else return 0.;
}
//_______________________________________________________________________________________
double GSystSet::Step(GSyst_t syst) const
{
  if ( this->IsIncluded(syst) ) return fSystematics.find(syst)->second->Step;
  else return 0.;
}
//_______________________________________________________________________________________
void GSystSet::SetCurValue(GSyst_t syst, double val)
{
  if ( this->IsIncluded(syst) ) {
    fSystematics[syst]->CurValue = val;
  }
}
//_______________________________________________________________________________________
void GSystSet::SetInitValue(GSyst_t syst, double val)
{
  if ( this->IsIncluded(syst) ) {
    fSystematics[syst]->InitValue = val;
  }
}
//_______________________________________________________________________________________
void GSystSet::SetRange(GSyst_t syst, double min, double max)
{
  if ( this->IsIncluded(syst) ) {
    fSystematics[syst]->MinValue = min;
    fSystematics[syst]->MaxValue = max;
  }
}
//_______________________________________________________________________________________
void GSystSet::SetStep(GSyst_t syst, double step)
{
  if ( this->IsIncluded(syst) ) {
    fSystematics[syst]->Step = step;
  }
}
//_______________________________________________________________________________________
void GSystSet::PrintSummary(void)
{
  LOG("ReW", pNOTICE) 
     << "Considering " << this->NIncluded() << " systematics";
				    
  vector<genie::rew::GSyst_t> svec = this->AllIncluded();

  unsigned int i=0;
  vector<genie::rew::GSyst_t>::const_iterator it = svec.begin();
  for( ; it != svec.end(); ++it) {
    GSyst_t syst = *it;
    LOG("ReW", pNOTICE) << "(" << i++ << ") : " << GSyst::AsString(syst);
  }
}
//_______________________________________________________________________________________
void GSystSet::Copy(const GSystSet & syst_set)
{
  return fSystematics.clear();

  map<GSyst_t, GSystInfo*>::const_iterator it = syst_set.fSystematics.begin();
  for( ; it != syst_set.fSystematics.end(); ++it) {
    GSyst_t     syst       = it->first;
    GSystInfo * syst_info  = it->second;

    double cur  = syst_info->CurValue;
    double init = syst_info->InitValue;
    double min  = syst_info->MinValue;
    double max  = syst_info->MaxValue;
    double step = syst_info->Step;

    this->Include      (syst);
    this->SetCurValue  (syst, cur);
    this->SetInitValue (syst, init);
    this->SetRange     (syst, min, max);
    this->SetStep      (syst, step);
  }
}
//_______________________________________________________________________________________
