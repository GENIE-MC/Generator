//_____________________________________________________________________________
/*!

\class    genie::nuvld::eDiffXSecTableFields

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004          
*/
//_____________________________________________________________________________

#include "DBUtils/eDiffXSecTableFields.h"

using namespace genie::nuvld;

ClassImp(eDiffXSecTableFields)

//____________________________________________________________________________
eDiffXSecTableFields::eDiffXSecTableFields() :
DBTableFields()
{
  AddField("name");
  AddField("measurement_tag");
  AddField("Sigma");
  AddField("Sigma_units");
  AddField("dSigma");
  AddField("E");
  AddField("E_units");
  AddField("EP");
  AddField("EP_units");
  AddField("Theta");
  AddField("Theta_units");
  AddField("Q2");
  AddField("Q2_units");
  AddField("W2");
  AddField("W2_units");
  AddField("Nu");
  AddField("Nu_units");
  AddField("Epsilon");
  AddField("Gamma");
  AddField("x");
}
//____________________________________________________________________________
eDiffXSecTableFields::eDiffXSecTableFields(const DBTableFields * fields):
DBTableFields(fields)
{

}
//____________________________________________________________________________
eDiffXSecTableFields::~eDiffXSecTableFields()
{

}
//____________________________________________________________________________


