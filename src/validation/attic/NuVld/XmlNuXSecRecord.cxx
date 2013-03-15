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

#include <cstdlib>
#include <string>
#include <typeinfo>

#include "ValidationTools/NuVld/XmlNuXSecRecord.h"

using std::string;
using std::endl;
using std::cerr;
using std::cout;

namespace genie {
namespace nuvld {
  
//____________________________________________________________________________
ostream & operator<<(ostream & stream, const XmlNuXSecRecord & rec_xsec) 
{  
  try {
   const XmlRecordStructure & rec_str = 
                            dynamic_cast<const XmlRecordStructure &> (rec_xsec);
   stream << rec_str;
  } 
  catch( std::bad_cast ) {
   stream 
        << "warning: could not dynamic_cast to a reference of the base object"
        << endl;
  }
   return stream;
}
//____________________________________________________________________________
XmlNuXSecRecord * XmlNuXSecRecord::_instance = 0; 
//____________________________________________________________________________
XmlNuXSecRecord::XmlNuXSecRecord() : XmlRecordStructure()
{
  string path = 
     string( getenv("GENIE") ) + 
     string("/src/ValidationTools/NuVld/rec/record_v_xsec.elements");

  ReadElements(path.c_str());
  
  _instance = 0;
}
//____________________________________________________________________________
XmlNuXSecRecord::~XmlNuXSecRecord()
{
  _instance = 0;
} 
//____________________________________________________________________________
XmlNuXSecRecord * XmlNuXSecRecord::Instance()
{
  if(_instance == 0) _instance = new XmlNuXSecRecord();
  return _instance;
}
//____________________________________________________________________________

} // nuvld namespace
} // genie namespace
