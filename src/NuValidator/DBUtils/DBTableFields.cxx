//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBTableFields

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004          
*/
//_____________________________________________________________________________

#include <cassert>

#include "DBUtils/DBTableFields.h"

using std::endl;
using namespace genie::nuvld;

ClassImp(DBTableFields)

//____________________________________________________________________________
namespace genie {
 namespace nuvld {
    ostream & operator << (ostream & stream, const DBTableFields & fields)
    {
      fields.Print(stream);
      return stream;
    }
 }
}    
//____________________________________________________________________________
DBTableFields::DBTableFields()
{

}
//____________________________________________________________________________
DBTableFields::DBTableFields(const DBTableFields * fields)
{
  vector<string>::const_iterator iter;

  for(iter = fields->_fields.begin(); 
                      iter != fields->_fields.end(); ++iter)
                                                      this->AddField( *iter );
}
//____________________________________________________________________________
DBTableFields::~DBTableFields()
{

}
//____________________________________________________________________________
unsigned int DBTableFields::NFields(void) const
{
  return _fields.size();
}
//____________________________________________________________________________
string DBTableFields::FieldName(unsigned int ifield) const
{
  return _fields[ifield];
}
//____________________________________________________________________________
unsigned int DBTableFields::FieldPos(string field_name) const
{  
  for(unsigned int i=0; i<_fields.size(); i++)
                          if( _fields[i].compare(field_name) == 0 ) return i;

  assert( false );
}
//____________________________________________________________________________
void DBTableFields::Print(ostream & stream) const
{
  vector<string>::const_iterator field_iter;

  for(field_iter = _fields.begin();
                          field_iter != _fields.end(); ++field_iter)
                                                stream << *field_iter << endl;
}
//____________________________________________________________________________
void DBTableFields::AddField(string field)
{
  _fields.push_back(field);
}
//____________________________________________________________________________
string DBTableFields::operator() (unsigned int ifield) const
{
  return FieldName(ifield);
}
//____________________________________________________________________________
