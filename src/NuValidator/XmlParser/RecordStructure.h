//_____________________________________________________________________________
/*!

\class    genie::nuvld::RecordStructure

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

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
  
class RecordStructure
{
public:

   friend ostream & operator << (ostream & stream, const RecordStructure & recbase);

   const map<string, vector<string> > & get_structure(void) const { return _elements; }

protected:

  RecordStructure();
  RecordStructure(const char * filename);
  virtual ~RecordStructure();

  int ReadElements(const char * filename);

  map<string, vector<string> > _elements;
};

} // nuvld namespace
} // genie namespace

#endif // _RECORD_STRUCTURE_H_
