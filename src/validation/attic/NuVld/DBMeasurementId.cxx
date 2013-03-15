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

#include <sstream>

#include <TSQLResult.h>
#include <TSQLRow.h>

#include  "ValidationTools/NuVld/DBMeasurementId.h"

using std::ostringstream;
using std::endl;

using namespace genie::nuvld;

//______________________________________________________________________________
namespace genie {
 namespace nuvld {
   ostream & operator << (ostream & stream, const DBMeasurementId & id)
   {
     id.Print(stream);  return stream;
   }
 }
}
//______________________________________________________________________________
DBMeasurementId::DBMeasurementId()
{
  _experiment      = "";
  _measurement_tag = "";
}
//______________________________________________________________________________
DBMeasurementId::DBMeasurementId(const DBMeasurementId * mid)
{
  _experiment      = mid->_experiment;
  _measurement_tag = mid->_measurement_tag;

  vector<string>::const_iterator str_iter;

  for(str_iter = mid->_author.begin();
                     str_iter != mid->_author.end(); ++str_iter)
                                                _author.push_back( *str_iter );  
  for(str_iter = mid->_journal.begin();
                     str_iter != mid->_journal.end(); ++str_iter)
                                               _journal.push_back( *str_iter );
  for(str_iter = mid->_year.begin();
                     str_iter != mid->_year.end(); ++str_iter)
                                                  _year.push_back( *str_iter );                                                
}
//______________________________________________________________________________
DBMeasurementId::DBMeasurementId(string experiment, string measurement_tag)
{
  _experiment      = experiment;
  _measurement_tag = measurement_tag;
}
//______________________________________________________________________________
DBMeasurementId::~DBMeasurementId()
{

}
//______________________________________________________________________________
string DBMeasurementId::Experiment(void) const
{
  return _experiment;
}
//______________________________________________________________________________
string DBMeasurementId::XmlMeasurementTag(void) const
{
  return _measurement_tag;
}
//______________________________________________________________________________
unsigned int DBMeasurementId::NRefs(void) const
{
  return _author.size();
}
//______________________________________________________________________________
string DBMeasurementId::Reference(unsigned int iref) const
{
  if( iref < NRefs() ) {

     ostringstream ref;
     ref << _experiment << ": "
                 << Author(iref) << "; " << Journal(iref) << ", " << Year(iref);
     return ref.str();

  } else
     return "----";
}
//______________________________________________________________________________
string DBMeasurementId::Author(unsigned int iref) const
{
  if( iref < NRefs() ) return _author[iref];
  else return "----";
}
//______________________________________________________________________________
string DBMeasurementId::Journal(unsigned int iref) const
{
  if( iref < NRefs() ) return _journal[iref];
  else return "----";
}
//______________________________________________________________________________
string DBMeasurementId::Year(unsigned int iref) const
{
  if( iref < NRefs() ) return _year[iref];
  else return "----";
}
//______________________________________________________________________________
void DBMeasurementId::Print(ostream & stream) const
{
  stream << "experiment: " << _experiment
                          << " / measurement_tag: " << _measurement_tag << endl;

  unsigned int nrefs = NRefs();

  for(unsigned int iref = 0; iref < nrefs; iref++) {
    
     stream << "REF/author   = " <<  Author(iref)    << endl;
     stream << "REF/journal  = " <<  Journal(iref)   << endl;
     stream << "REF/year     = " <<  Year(iref)      << endl;
     stream << "legend entry = " <<  Reference(iref) << endl;
  }
}
//______________________________________________________________________________
