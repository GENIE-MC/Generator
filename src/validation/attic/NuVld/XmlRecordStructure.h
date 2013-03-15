//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlRecordStructure

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RECORD_STRUCTURE_H_
#define _RECORD_STRUCTURE_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>

using std::ostream;
using std::string;
using std::vector;
using std::map;

namespace genie {
namespace nuvld {
  
class XmlRecordStructure
{
public:
   friend ostream & operator << (ostream & stream, const XmlRecordStructure & recbase);

   const map<string, vector<string> > & get_structure(void) const { return _elements; }

protected:
  XmlRecordStructure();
  XmlRecordStructure(const char * filename);
  virtual ~XmlRecordStructure();

  int ReadElements(const char * filename);

  map<string, vector<string> > _elements;
};

} // nuvld namespace
} // genie namespace

#endif // _RECORD_STRUCTURE_H_
