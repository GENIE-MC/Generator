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

#include "ValidationTools/NuVld/DBNuXSecTableRow.h"
#include "ValidationTools/NuVld/DBNuXSecTableFields.h"
#include "Messenger/Messenger.h"

using std::endl;
using namespace genie::nuvld;

ClassImp(DBNuXSecTableRow)

namespace genie {
namespace nuvld {

//____________________________________________________________________________
ostream & operator << (ostream & stream, const DBNuXSecTableRow & row)
{
  row.Print(stream);
  return stream;
}
//____________________________________________________________________________
DBNuXSecTableRow::DBNuXSecTableRow() :
DBTableRow()
{

}
//____________________________________________________________________________
DBNuXSecTableRow::DBNuXSecTableRow(TSQLRow * row) :
DBTableRow(new DBNuXSecTableFields, row)
{

}
//____________________________________________________________________________
DBNuXSecTableRow::DBNuXSecTableRow(const DBTableRow * db_row) :
DBTableRow(db_row)
{

}
//____________________________________________________________________________
DBNuXSecTableRow::~DBNuXSecTableRow()
{

}
//____________________________________________________________________________
string DBNuXSecTableRow::Experiment(void) const
{
  return this->Field("name");
}
//____________________________________________________________________________
string DBNuXSecTableRow::XmlMeasurementTag(void) const
{
  return this->Field("measurement_tag");
}
//____________________________________________________________________________
double DBNuXSecTableRow::XSec(void) const
{
  double xsec   = atof( this->Field("xsec").c_str() );
  double energy = E();
  string norm   = XSecNorm();
  string units  = XSecUnits();

  ApplyUnitsFactor (xsec, units);
  UndoXSecNorm     (xsec, energy, norm);

  LOG("NuVldDBI", pINFO) << "xsec = " << xsec << " (E = " << energy << ")";

  return xsec;
}
//____________________________________________________________________________
double DBNuXSecTableRow::StatErrP(void) const
{
  double xsec_err = atof( this->Field("stat_err_p").c_str() );

  double energy = this->E();
  string norm   = this->XSecNorm();
  string units  = this->XSecUnits();

  this->ApplyUnitsFactor (xsec_err, units);
  this->UndoXSecNorm     (xsec_err, energy, norm);

  // LOG("NuVldDBI", pINFO) 
  //      << "+dxsec[stat] = " << xsec_err << " (E = " << energy << ")";

  return xsec_err;
}
//____________________________________________________________________________
double DBNuXSecTableRow::StatErrM(void) const
{
  double xsec_err = atof( this->Field("stat_err_m").c_str() );

  double energy = this->E();
  string norm   = this->XSecNorm();
  string units  = this->XSecUnits();

  this->ApplyUnitsFactor (xsec_err, units);
  this->UndoXSecNorm     (xsec_err, energy, norm);

  //LOG("NuVldDBI", pINFO) 
  //      << "-dxsec[stat] = " << xsec_err << " (E = " << energy << ")";

  return xsec_err;
}
//____________________________________________________________________________
double DBNuXSecTableRow::SystErrP(void) const
{
  double xsec_err = atof( this->Field("syst_err_p").c_str() );

  if(this->IsPercentageErr(this->SystErrType()))
  {
    double xsec = this->XSec();
    xsec_err = xsec_err * xsec / 100.;
  } else {
    double energy = this->E();
    string norm   = this->XSecNorm();
    string units  = this->XSecUnits();

    this->ApplyUnitsFactor (xsec_err, units);
    this->UndoXSecNorm     (xsec_err, energy, norm);
  }

  //LOG("NuVldDBI", pINFO) 
  //      << "+dxsec[syst] = " << xsec_err << " (E = " << this->E() << ")";

  return xsec_err;
}
//____________________________________________________________________________
double DBNuXSecTableRow::SystErrM(void) const
{
  double xsec_err = atof( this->Field("syst_err_m").c_str() );

  if(this->IsPercentageErr(this->SystErrType()))
  {
    double xsec = this->XSec();
    xsec_err = xsec_err * xsec / 100.;
  } else {
    double energy = this->E();
    string norm   = this->XSecNorm();
    string units  = this->XSecUnits();

    this->ApplyUnitsFactor (xsec_err, units);
    this->UndoXSecNorm     (xsec_err, energy, norm);
  }

  //LOG("NuVldDBI", pINFO) 
  //      << "-dxsec[syst] = " << xsec_err << " (E = " << this->E() << ")";

  return xsec_err;
}
//____________________________________________________________________________
double DBNuXSecTableRow::ErrP(void) const
{
  double stat = this->StatErrP();
  double syst = this->SystErrP();

  return TMath::Sqrt( stat*stat + syst*syst );
}
//____________________________________________________________________________
double DBNuXSecTableRow::ErrM(void) const
{
  double stat = this->StatErrM();
  double syst = this->SystErrM();

  return TMath::Sqrt( stat*stat + syst*syst );
}
//____________________________________________________________________________
string DBNuXSecTableRow::XSecNorm(void) const
{
  return this->Field("xsec_norm");
}
//____________________________________________________________________________
string DBNuXSecTableRow::XSecUnits(void) const
{
  return this->Field("xsec_units");
}
//____________________________________________________________________________
string DBNuXSecTableRow::StatErrType(void) const
{
  return this->Field("stat_err_type");
}
//____________________________________________________________________________
string DBNuXSecTableRow::SystErrType(void) const
{
  return this->Field("syst_err_type");
}
//____________________________________________________________________________
double DBNuXSecTableRow::E(void) const
{
  return atof( this->Field("E").c_str() );
}
//____________________________________________________________________________
double DBNuXSecTableRow::Emin(void) const
{
  return atof( this->Field("E_min").c_str() );
}
//____________________________________________________________________________
double DBNuXSecTableRow::Emax(void) const
{
  return atof( this->Field("E_max").c_str() );
}
//____________________________________________________________________________
double DBNuXSecTableRow::dE(void) const
{
  return ( Emax() - Emin() );
}
//____________________________________________________________________________
void DBNuXSecTableRow::UndoXSecNorm(double & xsec, double E, string norm) const
{
// Data base entries are in xsec or xsec/E
// Compute xsec for xsec/E entries

  if (norm.compare("E") == 0)  xsec *= E;
}
//____________________________________________________________________________
void DBNuXSecTableRow::ApplyUnitsFactor(double & xsec, string factor) const
{
   static double cm2     = 1;
   static double GeV     = 1;
   static double nucleon = 1;

   double d = atof( &(factor.c_str())[1]);

   if      (factor.compare("*cm2") == 0)     xsec *= cm2;
   else if (factor.compare("/GeV") == 0)     xsec /= GeV;
   else if (factor.compare("/nucleon") == 0) xsec /= nucleon;
   else                                      xsec *= (1e+38 * d);
}
//____________________________________________________________________________
bool DBNuXSecTableRow::IsPercentageErr(string err_type) const
{
  LOG("NuVldDBI", pDEBUG) << err_type;

  bool isp = (err_type.find("%")    != string::npos) &&
             (err_type.find("*")    != string::npos) &&
             (err_type.find("xsec") != string::npos);

  if(isp) {
      LOG("NuVldDBI", pDEBUG) 
         << "Error is given as fraction of cross section value";
  }
  return isp;
}
//____________________________________________________________________________
void DBNuXSecTableRow::Print (ostream & stream) const
{
  stream << "E = " << this->Field("E") << " " << this->Field("E_units")
         << " --> xsec " << this->Field("xsec")
         << " " << this->Field("xsec_units")
         << endl;
}
//____________________________________________________________________________

} // nuvld namespace
} // genie namespace

