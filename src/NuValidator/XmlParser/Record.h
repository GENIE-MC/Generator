//_____________________________________________________________________________
/*!

\class    genie::nuvld::Record

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _RECORD_H_
#define _RECORD_H_

#include "RecordBase.h"

namespace genie {
namespace nuvld {
  
template<class T> class Record;
template<class T> ostream & operator << (ostream & stream, const Record<T> & rec);

template<class T> class Record : public RecordBase
{
public:
  Record();
  ~Record();

  void                 PrintStructure(void)          const;
  const vector<string> GetElements(void)             const;
  const vector<string> GetAttributes(string element) const;

  void         Add(string key, string value);
  const string Get(string key) const;

  friend ostream& operator<< <T>(ostream& stream, const Record<T> & rec);

private:

  T * _structure;
};

} // nuvld namespace
} // genie namespace

#endif // _RECORD_H_
