//____________________________________________________________________________
/*!

\class    genie::BaryonResList

\brief    Encapsulates a list of baryon resonances.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResList.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"

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

}
//____________________________________________________________________________
BaryonResList::~BaryonResList()
{

}
//____________________________________________________________________________
unsigned int BaryonResList::NResonances(void) const
{
  return fResVec->size();
}
//____________________________________________________________________________
string BaryonResList::ResonanceName(unsigned int ires) const
{
  if( ires >= 0 & ires < NResonances() ) {

    return res_utils::AsString( (*fResVec)[ires] );

  } else {

    SLOG("BaryonResList", pERROR) << "*** ires: " << ires
                           << " outside limits: [0, " << NResonances() << "]";
  }

  return "-";
}
//____________________________________________________________________________
Resonance_t BaryonResList::ResonanceId(unsigned int ires) const
{
  if( ires >= 0 & ires < NResonances() ) {

    return (*fResVec)[ires];

  } else {

    SLOG("BaryonResList", pERROR) << "*** ires: " << ires
                           << " outside limits: [0, " << NResonances() << "]";
  }

  return kNoResonance;
}
//____________________________________________________________________________
int BaryonResList::ResonancePdgCode(unsigned int ires) const
{
  return 0;
}
//____________________________________________________________________________
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

    Resonance_t res = res_utils::FromString( (*riter).c_str() );

    if( res == kNoResonance ) {

        SLOG("BaryonResList", pERROR) << "*** Unknown resonance: " << *riter;

    } else (*fResVec)[ires++] = res;
  }
}
//____________________________________________________________________________
void BaryonResList::Print(ostream & stream) const
{
  stream << "\n [-] Resonance List\n";

  vector<Resonance_t>::const_iterator riter;

  for(riter = fResVec->begin(); riter != fResVec->end(); ++riter) {

       stream << "  |--> RES: " << res_utils::AsString(*riter) << endl;
  }
}
//____________________________________________________________________________


