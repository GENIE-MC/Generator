//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Aug 01, 2003

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency.
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "ValidationTools/NuVld/XmlRecordBase.h"

using std::endl;

namespace genie {
namespace nuvld {
  
//_____________________________________________________________________________
ostream & operator << (ostream & stream, const XmlRecordBase & rec_base)
{
  map<string, string>::const_iterator it = rec_base._rec.begin();

  stream << endl;
  for( ; it != rec_base._rec.end(); ++it) {
     stream << it->first << " : " << it->second << endl;
  }
  return stream;
}
//_____________________________________________________________________________
XmlRecordBase::XmlRecordBase()
{

}
//_____________________________________________________________________________
XmlRecordBase::~XmlRecordBase()
{

}
//_____________________________________________________________________________

} // nuvld namespace
} // genie namespace
