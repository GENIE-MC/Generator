//_____________________________________________________________________________
/*!

\class    genie::nuvld::eDiffXSecTableRow

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#include "DBUtils/eDiffXSecTableRow.h"
#include "DBUtils/eDiffXSecTableFields.h"

using std::endl;
using namespace genie::nuvld;

ClassImp(eDiffXSecTableRow)

namespace genie {
namespace nuvld {

//____________________________________________________________________________
ostream & operator << (ostream & stream, const eDiffXSecTableRow & row)
{
  row.Print(stream);
  return stream;
}
//____________________________________________________________________________
eDiffXSecTableRow::eDiffXSecTableRow() :
DBTableRow()
{

}
//____________________________________________________________________________
eDiffXSecTableRow::eDiffXSecTableRow(TSQLRow * row) :
DBTableRow(new eDiffXSecTableFields, row)
{

}
//____________________________________________________________________________
eDiffXSecTableRow::eDiffXSecTableRow(const DBTableRow * db_row) :
DBTableRow(db_row)
{

}
//____________________________________________________________________________
eDiffXSecTableRow::~eDiffXSecTableRow()
{

}
//____________________________________________________________________________
string eDiffXSecTableRow::Experiment(void) const
{
  return Field("name");
}
//____________________________________________________________________________
string eDiffXSecTableRow::MeasurementTag(void) const
{
  return Field("measurement_tag");
}
//____________________________________________________________________________
string eDiffXSecTableRow::SigmaUnits(void) const
{
  return Field("Sigma_units");
}
//____________________________________________________________________________
string eDiffXSecTableRow::EUnits(void) const
{
  return Field("E_units");
}
//____________________________________________________________________________
string eDiffXSecTableRow::EPUnits(void) const
{
  return Field("EP_units");
}
//____________________________________________________________________________
string eDiffXSecTableRow::ThetaUnits(void) const
{
  return Field("Theta_units");
}
//____________________________________________________________________________
string eDiffXSecTableRow::Q2Units(void) const
{
  return Field("Q2_units");
}
//____________________________________________________________________________
string eDiffXSecTableRow::W2Units(void) const
{
  return Field("W2_units");
}
//____________________________________________________________________________
string eDiffXSecTableRow::NuUnits(void) const
{
  return Field("Nu_units");
}
//____________________________________________________________________________
double eDiffXSecTableRow::Sigma(void) const
{
  return atof( Field("Sigma").c_str() );
}
//____________________________________________________________________________
double eDiffXSecTableRow::dSigma(void) const
{
  return atof( Field("dSigma").c_str() );
}
//____________________________________________________________________________
double eDiffXSecTableRow::E(void) const
{
  return atof( Field("E").c_str() );
}
//____________________________________________________________________________
double eDiffXSecTableRow::EP(void) const
{
  return atof( Field("EP").c_str() );
}
//____________________________________________________________________________
double eDiffXSecTableRow::Theta(void) const
{
  return atof( Field("Theta").c_str() );
}
//____________________________________________________________________________
double eDiffXSecTableRow::Q2(void) const
{
  return atof( Field("Q2").c_str() );
}
//____________________________________________________________________________
double eDiffXSecTableRow::W2(void) const
{
  return atof( Field("W2").c_str() );
}
//____________________________________________________________________________
double eDiffXSecTableRow::Nu(void) const
{
  return atof( Field("Nu").c_str() );
}
//____________________________________________________________________________
double eDiffXSecTableRow::Epsilon(void) const
{
  return atof( Field("Epsilon").c_str() );
}
//____________________________________________________________________________
double eDiffXSecTableRow::Gamma(void) const
{
  return atof( Field("Gamma").c_str() );
}
//____________________________________________________________________________
double eDiffXSecTableRow::x(void) const
{
  return atof( Field("x").c_str() );
}
//____________________________________________________________________________
void eDiffXSecTableRow::Print(ostream & stream) const
{
  stream 
    << " E= "      << this->Field("E")     << " " << this->Field("E_units")
    << " EP= "     << this->Field("EP")    << " " << this->Field("EP_units")
    << " Theta= "  << this->Field("Theta") << " " << this->Field("Theta_units")         
    << " -> sig= " << this->Field("Sigma") << " " << this->Field("Sigma_units")
    << endl;  
}
//____________________________________________________________________________

} // nuvld namespace
} // genie namespace


