//_____________________________________________________________________________
/*!

\class    genie::nuvld::vXSecTableRow

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004          
*/
//_____________________________________________________________________________

#include <TMath.h>

#include "DBUtils/vXSecTableRow.h"
#include "DBUtils/vXSecTableFields.h"

using std::endl;
using namespace genie::nuvld;

ClassImp(vXSecTableRow)

namespace genie {
namespace nuvld {

//____________________________________________________________________________
ostream & operator << (ostream & stream, const vXSecTableRow & row)
{
  stream << "E = " << row.Field("E") << " " << row.Field("E_units")
         << " --> xsec " << row.Field("xsec") << " " << row.Field("xsec_units")
         << endl;
         
  return stream;
}
//____________________________________________________________________________
vXSecTableRow::vXSecTableRow() :
DBTableRow()
{

}
//____________________________________________________________________________
//vXSecTableRow::vXSecTableRow(const DBTableFields * fields, TSQLRow * row) :
//DBTableRow(fields, row)
//{
   
//}
vXSecTableRow::vXSecTableRow(TSQLRow * row) :
DBTableRow(new vXSecTableFields, row)
{

}
//____________________________________________________________________________
vXSecTableRow::vXSecTableRow(const DBTableRow * db_row) :
DBTableRow(db_row)
{

}
//____________________________________________________________________________
vXSecTableRow::~vXSecTableRow()
{

}
//____________________________________________________________________________
string vXSecTableRow::Experiment(void) const
{
  return Field("name");
}
//____________________________________________________________________________
string vXSecTableRow::MeasurementTag(void) const
{
  return Field("measurement_tag");
}
//____________________________________________________________________________
double vXSecTableRow::XSec(void) const
{
  double xsec   = atof( Field("xsec").c_str() );  
  double energy = E();  
  string norm   = XSecNorm();
  string units  = XSecUnits();
  
  ApplyUnitsFactor (xsec, units);   
  UndoXSecNorm     (xsec, energy, norm);

  return xsec;
}
//____________________________________________________________________________
double vXSecTableRow::StatErrP(void) const 
{
  double xsec_err = atof( Field("stat_err_p").c_str() );  

  double energy = E();
  string norm   = XSecNorm();
  string units  = XSecUnits();

  ApplyUnitsFactor (xsec_err, units);
  UndoXSecNorm     (xsec_err, energy, norm);

  return xsec_err;
}
//____________________________________________________________________________
double vXSecTableRow::StatErrM(void) const
{
  double xsec_err = atof( Field("stat_err_m").c_str() );  

  double energy = E();
  string norm   = XSecNorm();
  string units  = XSecUnits();

  ApplyUnitsFactor (xsec_err, units);
  UndoXSecNorm     (xsec_err, energy, norm);

  return xsec_err;  
}
//____________________________________________________________________________
double vXSecTableRow::SystErrP(void) const
{
  return 0;
}
//____________________________________________________________________________
double vXSecTableRow::SystErrM(void) const
{
  return 0;
}
//____________________________________________________________________________
double vXSecTableRow::ErrP(void) const
{
  double stat = StatErrP();
  double syst = SystErrP();

  return TMath::Sqrt( stat*stat + syst*syst );
}
//____________________________________________________________________________
double vXSecTableRow::ErrM(void) const
{
  double stat = StatErrM();
  double syst = SystErrM();

  return TMath::Sqrt( stat*stat + syst*syst );
}
//____________________________________________________________________________
string vXSecTableRow::XSecNorm(void) const
{
  return Field("xsec_norm");
}
//____________________________________________________________________________
string vXSecTableRow::XSecUnits(void) const
{
  return Field("xsec_units");
}
//____________________________________________________________________________
double vXSecTableRow::E(void) const
{
  return atof( Field("E").c_str() );
}
//____________________________________________________________________________
double vXSecTableRow::Emin(void) const
{
  return atof( Field("E_min").c_str() );
}
//____________________________________________________________________________
double vXSecTableRow::Emax(void) const
{
  return atof( Field("E_max").c_str() );
}
//____________________________________________________________________________
double vXSecTableRow::dE(void) const
{
  return ( Emax() - Emin() );
}

//____________________________________________________________________________
void vXSecTableRow::UndoXSecNorm(double & xsec, double E, string norm) const
{
// Data base entries are in xsec or xsec/E
// Compute xsec for xsec/E entries

  if (norm.compare("E") == 0)  xsec *= E;
}
//____________________________________________________________________________
void vXSecTableRow::ApplyUnitsFactor(double & xsec, string factor) const
{
   static double cm2     = 1;
   static double GeV     = 1;
   static double nucleon = 1;

   if      (factor.compare("*cm2") == 0)     xsec *= cm2;
   else if (factor.compare("/GeV") == 0)     xsec /= GeV;
   else if (factor.compare("/nucleon") == 0) xsec /= nucleon;
   else                                      xsec *= (1e+38 * atof( &(factor.c_str())[1]));
}
//____________________________________________________________________________

} // nuvld namespace
} // genie namespace

