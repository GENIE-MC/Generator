//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlRecord

\brief                                                           

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RECORD_H_
#define _RECORD_H_

#include "ValidationTools/NuVld/XmlRecordBase.h"

namespace genie {
namespace nuvld {
  
template<class T> class XmlRecord;
template<class T> ostream & operator << (ostream & stream, const XmlRecord<T> & rec);

template<class T> class XmlRecord : public XmlRecordBase
{
public:
  XmlRecord();
  ~XmlRecord();

  void                 PrintStructure(void)          const;
  const vector<string> GetElements(void)             const;
  const vector<string> GetAttributes(string element) const;

  void         Add(string key, string value);
  const string Get(string key) const;

  friend ostream& operator<< <T>(ostream& stream, const XmlRecord<T> & rec);

private:

  T * _structure;
};

} // nuvld namespace
} // genie namespace

#endif // _RECORD_H_
