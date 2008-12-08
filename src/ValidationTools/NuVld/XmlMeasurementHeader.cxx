//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - Aug 01, 2003

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency.
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "ValidationTools/NuVld/XmlMeasurementHeader.h"

namespace genie {
namespace nuvld {

//___________________________________________________________________________
ostream & operator << (ostream & stream, const XmlMeasurementHeader & header)
{
  stream << endl;
  stream << "---------------------- printing measurement header: " << endl;

  const vector<XmlCitation *> & refs = header.GetRefs();

  vector<XmlCitation *>::const_iterator ref_iter;
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
  stream << " |-> error status.....: " << header._err_status     << endl;
  stream << " |-> comment..........: " << header._comment        << endl;

  return stream;
}
//___________________________________________________________________________
XmlMeasurementHeader::XmlMeasurementHeader() :
_observable(""),
_target(""),
_reaction(""),
_A("1"),
_exposure(""),
_exposure_units(""),
_data_source(""),
_tag(""),
_npoints(""),
_err_status("1"),
_comment("")
{
  _refs = new vector<XmlCitation *>;
}
//___________________________________________________________________________
XmlMeasurementHeader::XmlMeasurementHeader(const XmlMeasurementHeader & /*header*/)
{

}
//___________________________________________________________________________
void XmlMeasurementHeader::Add(string key, string value)
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
  else if (key.compare("error_status")   == 0) _err_status     = value;
  else if (key.compare("comment")        == 0) _comment        = value;
  else
      cerr << "XmlMeasurementHeader::add(string, string): unknown key" << endl;
}
//___________________________________________________________________________
void XmlMeasurementHeader::Add(XmlCitation * ref)
{
  _refs->push_back(ref);
}
//___________________________________________________________________________

} // nuvld namespace
} // genie namespace
