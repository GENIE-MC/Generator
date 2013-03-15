//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include <TMath.h>

#include "ValidationTools/NuVld/DBSFTableRow.h"
#include "ValidationTools/NuVld/DBSFTableFields.h"

using std::endl;
using namespace genie::nuvld;

ClassImp(DBSFTableRow)

namespace genie {
namespace nuvld {

//____________________________________________________________________________
ostream & operator << (ostream & stream, const DBSFTableRow & row)
{
  row.Print(stream);         
  return stream;
}
//____________________________________________________________________________
DBSFTableRow::DBSFTableRow() :
DBTableRow()
{

}
//____________________________________________________________________________
DBSFTableRow::DBSFTableRow(TSQLRow * row) :
DBTableRow(new DBSFTableFields, row)
{

}
//____________________________________________________________________________
DBSFTableRow::DBSFTableRow(const DBTableRow * db_row) :
DBTableRow(db_row)
{

}
//____________________________________________________________________________
DBSFTableRow::~DBSFTableRow()
{

}
//____________________________________________________________________________
string DBSFTableRow::Experiment(void) const
{
  return Field("name");
}
//____________________________________________________________________________
string DBSFTableRow::XmlMeasurementTag(void) const
{
  return Field("measurement_tag");
}
//____________________________________________________________________________
double DBSFTableRow::SF(void) const
{
  return atof( Field("sf").c_str() );
}
//____________________________________________________________________________
double DBSFTableRow::StatErrP(void) const 
{
  return atof( Field("stat_err_p").c_str() );  
}
//____________________________________________________________________________
double DBSFTableRow::StatErrM(void) const
{
  return atof( Field("stat_err_m").c_str() );  
}
//____________________________________________________________________________
double DBSFTableRow::SystErrP(void) const
{
  return atof( Field("syst_err_p").c_str() );
}
//____________________________________________________________________________
double DBSFTableRow::SystErrM(void) const
{
  return atof( Field("syst_err_m").c_str() );
}
//____________________________________________________________________________
double DBSFTableRow::ErrP(void) const
{
  double stat = StatErrP();
  double syst = SystErrP();

  return TMath::Sqrt( stat*stat + syst*syst );
}
//____________________________________________________________________________
double DBSFTableRow::ErrM(void) const
{
  double stat = StatErrM();
  double syst = SystErrM();

  return TMath::Sqrt( stat*stat + syst*syst );
}
//____________________________________________________________________________
string DBSFTableRow::R(void) const
{
  return Field("R");
}
//____________________________________________________________________________
string DBSFTableRow::p(void) const
{
  return Field("p");
}
//____________________________________________________________________________
double DBSFTableRow::Q2(void) const
{
  return atof( Field("Q2").c_str() );
}
//____________________________________________________________________________
double DBSFTableRow::x(void) const
{
  return atof( Field("x").c_str() );
}
//____________________________________________________________________________
void DBSFTableRow::Print(ostream & stream) const
{
  stream << "x = " << this->x() 
         << ", Q2 = " << this->Q2() << ", SF = " << this->SF() << endl;
}
//____________________________________________________________________________

} // nuvld namespace
} // genie namespace

