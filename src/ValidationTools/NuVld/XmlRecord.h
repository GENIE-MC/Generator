//_____________________________________________________________________________
/*!

\class    genie::nuvld::XmlRecord

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

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
