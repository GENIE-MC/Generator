//_____________________________________________________________________________
/*!

\class    genie::nuvld::MeasurementId

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004          
*/
//_____________________________________________________________________________

#include <sstream>

#include <TSQLResult.h>
#include <TSQLRow.h>

#include "DBUtils/MeasurementId.h"

using std::ostringstream;
using std::endl;

using namespace genie::nuvld;

//______________________________________________________________________________
namespace genie {
 namespace nuvld {
   ostream & operator << (ostream & stream, const MeasurementId & id)
   {
     id.Print(stream);  return stream;
   }
 }
}
//______________________________________________________________________________
MeasurementId::MeasurementId()
{
  _experiment      = "";
  _measurement_tag = "";
}
//______________________________________________________________________________
MeasurementId::MeasurementId(const MeasurementId * mid)
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
MeasurementId::MeasurementId(string experiment, string measurement_tag)
{
  _experiment      = experiment;
  _measurement_tag = measurement_tag;
}
//______________________________________________________________________________
MeasurementId::~MeasurementId()
{

}
//______________________________________________________________________________
string MeasurementId::Experiment(void) const
{
  return _experiment;
}
//______________________________________________________________________________
string MeasurementId::MeasurementTag(void) const
{
  return _measurement_tag;
}
//______________________________________________________________________________
unsigned int MeasurementId::NRefs(void) const
{
  return _author.size();
}
//______________________________________________________________________________
string MeasurementId::Reference(unsigned int iref) const
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
string MeasurementId::Author(unsigned int iref) const
{
  if( iref < NRefs() ) return _author[iref];
  else return "----";
}
//______________________________________________________________________________
string MeasurementId::Journal(unsigned int iref) const
{
  if( iref < NRefs() ) return _journal[iref];
  else return "----";
}
//______________________________________________________________________________
string MeasurementId::Year(unsigned int iref) const
{
  if( iref < NRefs() ) return _year[iref];
  else return "----";
}
//______________________________________________________________________________
void MeasurementId::Print(ostream & stream) const
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
