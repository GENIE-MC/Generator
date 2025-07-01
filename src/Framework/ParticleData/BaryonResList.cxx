//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <algorithm>

#include "Framework/ParticleData//BaryonResList.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/StringUtils.h"

using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie
{
  ostream & operator<<(ostream & stream, const BaryonResList & res_list)
  {
     res_list.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
BaryonResList::BaryonResList()
{
  fResVec = 0;
}
//____________________________________________________________________________
BaryonResList::BaryonResList(const BaryonResList & res_list)
{
  fResVec = 0;
  this->Copy(res_list);
}
//____________________________________________________________________________
BaryonResList::~BaryonResList()
{
  if(fResVec) delete fResVec;
}
//____________________________________________________________________________
unsigned int BaryonResList::NResonances(void) const
{
  if(!fResVec) {
    SLOG("BaryonResList", pERROR) << "Null Resonance List";
    return 0;
  }
  return fResVec->size();
}
//____________________________________________________________________________
string BaryonResList::ResonanceName(unsigned int ires) const
{
  if(!fResVec) {
    SLOG("BaryonResList", pERROR) << "Null Resonance List";
    return "-";
  }
  if(ires >= this->NResonances() ) {
    SLOG("BaryonResList", pERROR) << "Resonance idx: " << ires
                   << " outside limits: [0, " << this->NResonances() << "]";
    return "-";
  }
  return utils::res::AsString( (*fResVec)[ires] );
}
//____________________________________________________________________________
Resonance_t BaryonResList::ResonanceId(unsigned int ires) const
{
  if(!fResVec) {
    SLOG("BaryonResList", pERROR) << "Null Resonance List";
    return kNoResonance;
  }
  if(ires >= this->NResonances() ) {
    SLOG("BaryonResList", pERROR) << "Resonance idx: " << ires
                   << " outside limits: [0, " << this->NResonances() << "]";
    return kNoResonance;
  }
  return (*fResVec)[ires];
}
//____________________________________________________________________________
int BaryonResList::ResonancePdgCode(unsigned int /*ires*/) const
{
  return 0;
}
//____________________________________________________________________________
bool BaryonResList::Find(Resonance_t res) const
{
  if(!fResVec) {
    SLOG("BaryonResList", pWARN) << "NULL resonance list!";
    return false;
  }
  int n = count(fResVec->begin(), fResVec->end(), res);
  if(n!=0) return true;
  return false;
}
//___________________________________________________________________________
void BaryonResList::DecodeFromNameList(string input_list, string delimiter)
{
  //-- remove all spaces in the input string coming from the XML config file

  string list = utils::str::FilterString(" ", input_list);

  vector<string> resonances = utils::str::Split(list, delimiter);

  SLOG("BaryonResList", pINFO) << list;
  SLOG("BaryonResList", pINFO) << resonances.size();

  if(fResVec) delete fResVec;
  fResVec = new vector<Resonance_t> (resonances.size());

  unsigned int ires = 0;
  vector<string>::const_iterator riter;
  for(riter = resonances.begin(); riter != resonances.end(); ++riter) {

    Resonance_t res = utils::res::FromString( (*riter).c_str() );
    if( res == kNoResonance ) {
        SLOG("BaryonResList", pERROR) << "*** Unknown resonance: " << *riter;
    } else (*fResVec)[ires++] = res;
  }
}
//____________________________________________________________________________
void BaryonResList::Clear(void)
{
  if(fResVec) fResVec->clear();
}
//____________________________________________________________________________
void BaryonResList::Copy(const BaryonResList & res_list)
{
  if(fResVec) fResVec->clear();

  unsigned int nres = res_list.NResonances();
  if(nres==0) return;

  if(!fResVec) fResVec = new vector<Resonance_t> (nres);

  for(unsigned int ires = 0; ires < nres; ires++) {
     (*fResVec)[ires] = res_list.ResonanceId(ires);
  }
}
//____________________________________________________________________________
void BaryonResList::Print(ostream & stream) const
{
  stream << "\n [-] Resonance List\n";

  vector<Resonance_t>::const_iterator riter;
  for(riter = fResVec->begin(); riter != fResVec->end(); ++riter) {
       stream << "  |--> RES: " << utils::res::AsString(*riter) << endl;
  }
}
//____________________________________________________________________________
auto BaryonResList::begin() noexcept -> typename vector<Resonance_t>::iterator
{
  return fResVec->begin();
}
//____________________________________________________________________________
auto BaryonResList::end() noexcept -> typename vector<Resonance_t>::iterator
{
  return fResVec->end();
}
//____________________________________________________________________________
auto BaryonResList::begin() const noexcept -> typename vector<Resonance_t>::const_iterator
{
  return fResVec->begin();
}
//____________________________________________________________________________
auto BaryonResList::end() const noexcept -> typename vector<Resonance_t>::const_iterator
{
  return fResVec->end();
}
//____________________________________________________________________________
auto BaryonResList::cbegin() const noexcept -> typename vector<Resonance_t>::const_iterator
{
  return fResVec->cbegin();
}
//____________________________________________________________________________
auto BaryonResList::cend() const noexcept -> typename vector<Resonance_t>::const_iterator
{
  return fResVec->cend();
}
//____________________________________________________________________________

