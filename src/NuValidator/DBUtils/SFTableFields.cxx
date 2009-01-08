//_____________________________________________________________________________
/*!

\class    genie::nuvld::SFTableFields

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004          
*/
//_____________________________________________________________________________

#include "DBUtils/SFTableFields.h"

using namespace genie::nuvld;

ClassImp(SFTableFields)

//____________________________________________________________________________
SFTableFields::SFTableFields() :
DBTableFields()
{
  AddField("name");
  AddField("measurement_tag");
  AddField("sf");
  AddField("stat_err_p");
  AddField("stat_err_m");
  AddField("syst_err_p");
  AddField("syst_err_m");
  AddField("R");
  AddField("p");
  AddField("x");
  AddField("Q2");
}
//____________________________________________________________________________
SFTableFields::SFTableFields(const DBTableFields * fields):
DBTableFields(fields)
{

}
//____________________________________________________________________________
SFTableFields::~SFTableFields()
{

}
//____________________________________________________________________________


