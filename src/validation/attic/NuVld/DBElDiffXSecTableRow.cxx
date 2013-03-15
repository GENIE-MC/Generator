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

#include <cstdlib>

#include "ValidationTools/NuVld/DBElDiffXSecTableRow.h"
#include "ValidationTools/NuVld/DBElDiffXSecTableFields.h"

using std::endl;
using namespace genie::nuvld;

ClassImp(DBElDiffXSecTableRow)

namespace genie {
namespace nuvld {

//____________________________________________________________________________
ostream & operator << (ostream & stream, const DBElDiffXSecTableRow & row)
{
  row.Print(stream);
  return stream;
}
//____________________________________________________________________________
DBElDiffXSecTableRow::DBElDiffXSecTableRow():
DBTableRow()
{

}
//____________________________________________________________________________
DBElDiffXSecTableRow::DBElDiffXSecTableRow(TSQLRow * row) :
DBTableRow(new DBElDiffXSecTableFields, row)
{

}
//____________________________________________________________________________
DBElDiffXSecTableRow::DBElDiffXSecTableRow(const DBTableRow * db_row) :
DBTableRow(db_row)
{

}
//____________________________________________________________________________
DBElDiffXSecTableRow::~DBElDiffXSecTableRow()
{

}
//____________________________________________________________________________
string DBElDiffXSecTableRow::Experiment(void) const
{
  return Field("name");
}
//____________________________________________________________________________
string DBElDiffXSecTableRow::XmlMeasurementTag(void) const
{
  return Field("measurement_tag");
}
//____________________________________________________________________________
string DBElDiffXSecTableRow::SigmaUnits(void) const
{
  return Field("Sigma_units");
}
//____________________________________________________________________________
string DBElDiffXSecTableRow::EUnits(void) const
{
  return Field("E_units");
}
//____________________________________________________________________________
string DBElDiffXSecTableRow::EPUnits(void) const
{
  return Field("EP_units");
}
//____________________________________________________________________________
string DBElDiffXSecTableRow::ThetaUnits(void) const
{
  return Field("Theta_units");
}
//____________________________________________________________________________
string DBElDiffXSecTableRow::Q2Units(void) const
{
  return Field("Q2_units");
}
//____________________________________________________________________________
string DBElDiffXSecTableRow::W2Units(void) const
{
  return Field("W2_units");
}
//____________________________________________________________________________
string DBElDiffXSecTableRow::NuUnits(void) const
{
  return Field("Nu_units");
}
//____________________________________________________________________________
double DBElDiffXSecTableRow::Sigma(void) const
{
  return atof( Field("Sigma").c_str() );
}
//____________________________________________________________________________
double DBElDiffXSecTableRow::dSigma(void) const
{
  return atof( Field("dSigma").c_str() );
}
//____________________________________________________________________________
double DBElDiffXSecTableRow::E(void) const
{
  return atof( Field("E").c_str() );
}
//____________________________________________________________________________
double DBElDiffXSecTableRow::EP(void) const
{
  return atof( Field("EP").c_str() );
}
//____________________________________________________________________________
double DBElDiffXSecTableRow::Theta(void) const
{
  return atof( Field("Theta").c_str() );
}
//____________________________________________________________________________
double DBElDiffXSecTableRow::Q2(void) const
{
  return atof( Field("Q2").c_str() );
}
//____________________________________________________________________________
double DBElDiffXSecTableRow::W2(void) const
{
  return atof( Field("W2").c_str() );
}
//____________________________________________________________________________
double DBElDiffXSecTableRow::Nu(void) const
{
  return atof( Field("Nu").c_str() );
}
//____________________________________________________________________________
double DBElDiffXSecTableRow::Epsilon(void) const
{
  return atof( Field("Epsilon").c_str() );
}
//____________________________________________________________________________
double DBElDiffXSecTableRow::Gamma(void) const
{
  return atof( Field("Gamma").c_str() );
}
//____________________________________________________________________________
double DBElDiffXSecTableRow::x(void) const
{
  return atof( Field("x").c_str() );
}
//____________________________________________________________________________
void DBElDiffXSecTableRow::Print(ostream & stream) const
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


