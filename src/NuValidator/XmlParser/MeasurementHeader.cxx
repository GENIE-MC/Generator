//_____________________________________________________________________________
/*!

\class    genie::nuvld::MeasurementHeader

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include "MeasurementHeader.h"

namespace genie {
namespace nuvld {
  
//___________________________________________________________________________
ostream & operator << (ostream & stream, const MeasurementHeader & header)
{
  stream << endl;
  stream << "---------------------- printing measurement header: " << endl;

  const vector<Citation *> & refs = header.GetRefs();
  
  vector<Citation *>::const_iterator ref_iter;
  for(ref_iter = refs.begin();
                ref_iter != refs.end(); ++ref_iter) stream << *(*ref_iter);

  stream << "[-] " << endl;
  stream << " |-> observable.......: " << header._observable     << endl;
  stream << " |-> target...........: " << header._target         << endl;
  stream << " |-> reaction.........: " << header._reaction       << endl;
  stream << " |-> A................: " << header._A              << endl;
  stream << " |-> exposure.........: " << header._exposure       << endl;
  stream << " |-> exposure units...: " << header._exposure_units << endl;
  stream << " |-> data source......: " << header._data_source    << endl;
  stream << " |-> tag..............: " << header._tag            << endl;
  stream << " |-> npoints..........: " << header._npoints        << endl;
  stream << " |-> comment..........: " << header._comment        << endl;

  return stream;
}
//___________________________________________________________________________
MeasurementHeader::MeasurementHeader() :
_observable(""),
_target(""),
_reaction(""),
_A("1"),
_exposure(""),
_exposure_units(""),
_data_source(""),
_tag(""),
_npoints(""),
_comment("")
{
  _refs = new vector<Citation *>;
}
//___________________________________________________________________________
MeasurementHeader::MeasurementHeader(const MeasurementHeader & header)
{

}
//___________________________________________________________________________
void MeasurementHeader::Add(string key, string value)
{
  if      (key.compare("observable")     == 0) _observable     = value;
  else if (key.compare("target")         == 0) _target         = value;
  else if (key.compare("reaction")       == 0) _reaction       = value;
  else if (key.compare("A")              == 0) _A              = value;
  else if (key.compare("exposure")       == 0) _exposure       = value;
  else if (key.compare("exposure_units") == 0) _exposure_units = value;
  else if (key.compare("data_source")    == 0) _data_source    = value;
  else if (key.compare("tag")            == 0) _tag            = value;
  else if (key.compare("npoints")        == 0) _npoints        = value;
  else if (key.compare("comment")        == 0) _comment        = value;
  else 
      cerr << "MeasurementHeader::add(string, string): unknown key" << endl;
}
//___________________________________________________________________________
void MeasurementHeader::Add(Citation * ref)
{
  _refs->push_back(ref);
}
//___________________________________________________________________________

} // nuvld namespace
} // genie namespace
