//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency.
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "Messenger/Messenger.h"
#include "ValidationTools/NuVld/XmlMeasurementHeader.h"

namespace genie {
namespace nuvld {

//___________________________________________________________________________
ostream & operator << (ostream & stream, const XmlMeasurementHeader & header)
{
  stream << endl;

  const vector<XmlCitation *> & refs = header.GetRefs();
  vector<XmlCitation *>::const_iterator ref_iter = refs.begin();
  for( ; ref_iter != refs.end(); ++ref_iter) {
     stream << *(*ref_iter);
  }

  stream << "[-] " << endl;
  stream << " |-> tag..............: " << header.fTag            << endl;
  stream << " |-> observable.......: " << header.fObservable     << endl;
  stream << " |-> target...........: " << header.fTarget         << endl;
  stream << " |-> reaction.........: " << header.fReaction       << endl;
  stream << " |-> A................: " << header.fA              << endl;
  stream << " |-> exposure.........: " << header.fExposure       << endl;
  stream << " |-> exposure units...: " << header.fExposureUnits  << endl;
  stream << " |-> data source......: " << header.fDataSource     << endl;
  stream << " |-> npoints..........: " << header.fNpoints        << endl;
  stream << " |-> error status.....: " << header.fErrStatus      << endl;
  stream << " |-> comment..........: " << header.fComment        << endl;

  return stream;
}
//___________________________________________________________________________
XmlMeasurementHeader::XmlMeasurementHeader() :
fTag           (""),
fObservable    (""),
fTarget        (""),
fReaction      (""),
fA             ("1"),
fExposure      (""),
fExposureUnits (""),
fDataSource    (""),
fNpoints       (""),
fErrStatus     ("1"),
fComment       ("")
{
  fRefs = new vector<XmlCitation *>;
}
//___________________________________________________________________________
XmlMeasurementHeader::XmlMeasurementHeader(const XmlMeasurementHeader& /*h*/)
{

}
//___________________________________________________________________________
void XmlMeasurementHeader::Add(string key, string value)
{
  if      (key.compare("tag")            == 0) fTag            = value;
  else if (key.compare("observable")     == 0) fObservable     = value;
  else if (key.compare("target")         == 0) fTarget         = value;
  else if (key.compare("reaction")       == 0) fReaction       = value;
  else if (key.compare("A")              == 0) fA              = value;
  else if (key.compare("exposure")       == 0) fExposure       = value;
  else if (key.compare("exposure_units") == 0) fExposureUnits  = value;
  else if (key.compare("data_source")    == 0) fDataSource     = value;
  else if (key.compare("npoints")        == 0) fNpoints        = value;
  else if (key.compare("error_status")   == 0) fErrStatus      = value;
  else if (key.compare("comment")        == 0) fComment        = value;
  else {
   LOG("NuVld", pERROR) << "Unknown key: " << key;
  }
}
//___________________________________________________________________________
void XmlMeasurementHeader::Add(XmlCitation * ref)
{
  fRefs->push_back(ref);
}
//___________________________________________________________________________

} // nuvld namespace
} // genie namespace
