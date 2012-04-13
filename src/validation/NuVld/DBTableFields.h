//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBTableFields

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _DB_TABLE_FIELDS_H_
#define _DB_TABLE_FIELDS_H_

#include <vector>
#include <string>
#include <ostream>

#include <TObject.h>

using std::vector;
using std::string;
using std::ostream;

namespace genie {
namespace nuvld {

class DBTableFields : public TObject
{
public:

  DBTableFields();
  DBTableFields(const DBTableFields * fields);
  ~DBTableFields();

  unsigned int NFields   (void)                const;  
  string       FieldName (unsigned int ifield) const;
  unsigned int FieldPos  (string field_name)   const;

  void Print(ostream & stream) const;

  string operator() (unsigned int ifield) const;
  
  friend ostream & operator << (ostream & stream, const DBTableFields & fields);
    
protected:  

  void AddField(string field);
  
  vector<string> _fields;

ClassDef(DBTableFields, 1)
};

} // nuvld namespace
} // genie namespace

#endif
