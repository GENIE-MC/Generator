//_____________________________________________________________________________
/*!

\class    genie::nuvld::vXSecTableFields

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#include "DBUtils/vXSecTableFields.h"

using namespace genie::nuvld;

ClassImp(vXSecTableFields)

//____________________________________________________________________________
vXSecTableFields::vXSecTableFields() :
DBTableFields()
{
  AddField("name");
  AddField("measurement_tag");
  AddField("xsec");
  AddField("stat_err_p");
  AddField("stat_err_m");
  AddField("syst_err_p");
  AddField("syst_err_m");
  AddField("xsec_units");
  AddField("xsec_norm");
  AddField("stat_err_type");
  AddField("syst_err_type");
  AddField("E");
  AddField("E_min");
  AddField("E_max");
  AddField("E_units");
  AddField("E_frame");
}
//____________________________________________________________________________
vXSecTableFields::vXSecTableFields(const DBTableFields * fields):
DBTableFields(fields)
{

}
//____________________________________________________________________________
vXSecTableFields::~vXSecTableFields()
{

}
//____________________________________________________________________________


