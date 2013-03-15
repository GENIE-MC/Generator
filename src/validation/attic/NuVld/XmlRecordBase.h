//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlRecordBase

\brief                                                           

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RECORD_BASE_H_
#define _RECORD_BASE_H_

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
  
class XmlRecordBase
{
public:
   XmlRecordBase();
   virtual ~XmlRecordBase();

   friend ostream & operator << (ostream & stream, const XmlRecordBase & rec_base);

   virtual void                 PrintStructure(void) const          = 0;
   virtual const vector<string> GetElements(void) const             = 0;
   virtual const vector<string> GetAttributes(string element) const = 0;

   virtual void                 Add(string key, string value)        = 0;
   virtual const string         Get(string key) const                = 0;   

protected:

  map<string, string > _rec;
};

} // nuvld namespace
} // genie namespace

#endif // _RECORD_BASE_H_
