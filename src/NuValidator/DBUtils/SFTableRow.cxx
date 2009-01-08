//_____________________________________________________________________________
/*!

\class    genie::nuvld::SFTableRow

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004          
*/
//_____________________________________________________________________________

#include <TMath.h>

#include "DBUtils/SFTableRow.h"
#include "DBUtils/SFTableFields.h"

using std::endl;
using namespace genie::nuvld;

ClassImp(SFTableRow)

namespace genie {
namespace nuvld {

//____________________________________________________________________________
ostream & operator << (ostream & stream, const SFTableRow & row)
{
  row.Print(stream);         
  return stream;
}
//____________________________________________________________________________
SFTableRow::SFTableRow() :
DBTableRow()
{

}
//____________________________________________________________________________
SFTableRow::SFTableRow(TSQLRow * row) :
DBTableRow(new SFTableFields, row)
{

}
//____________________________________________________________________________
SFTableRow::SFTableRow(const DBTableRow * db_row) :
DBTableRow(db_row)
{

}
//____________________________________________________________________________
SFTableRow::~SFTableRow()
{

}
//____________________________________________________________________________
string SFTableRow::Experiment(void) const
{
  return Field("name");
}
//____________________________________________________________________________
string SFTableRow::MeasurementTag(void) const
{
  return Field("measurement_tag");
}
//____________________________________________________________________________
double SFTableRow::SF(void) const
{
  return atof( Field("sf").c_str() );
}
//____________________________________________________________________________
double SFTableRow::StatErrP(void) const 
{
  return atof( Field("stat_err_p").c_str() );  
}
//____________________________________________________________________________
double SFTableRow::StatErrM(void) const
{
  return atof( Field("stat_err_m").c_str() );  
}
//____________________________________________________________________________
double SFTableRow::SystErrP(void) const
{
  return atof( Field("syst_err_p").c_str() );
}
//____________________________________________________________________________
double SFTableRow::SystErrM(void) const
{
  return atof( Field("syst_err_m").c_str() );
}
//____________________________________________________________________________
double SFTableRow::ErrP(void) const
{
  double stat = StatErrP();
  double syst = SystErrP();

  return TMath::Sqrt( stat*stat + syst*syst );
}
//____________________________________________________________________________
double SFTableRow::ErrM(void) const
{
  double stat = StatErrM();
  double syst = SystErrM();

  return TMath::Sqrt( stat*stat + syst*syst );
}
//____________________________________________________________________________
string SFTableRow::R(void) const
{
  return Field("R");
}
//____________________________________________________________________________
string SFTableRow::p(void) const
{
  return Field("p");
}
//____________________________________________________________________________
double SFTableRow::Q2(void) const
{
  return atof( Field("Q2").c_str() );
}
//____________________________________________________________________________
double SFTableRow::x(void) const
{
  return atof( Field("x").c_str() );
}
//____________________________________________________________________________
void SFTableRow::Print(ostream & stream) const
{
  stream << "x = " << this->x() 
         << ", Q2 = " << this->Q2() << ", SF = " << this->SF() << endl;
}
//____________________________________________________________________________

} // nuvld namespace
} // genie namespace

