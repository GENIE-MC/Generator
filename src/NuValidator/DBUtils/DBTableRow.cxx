//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBTableRow

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004          
*/
//_____________________________________________________________________________

#include "DBUtils/DBTableRow.h"

using std::endl;
using namespace genie::nuvld;

ClassImp(DBTableRow)

//____________________________________________________________________________
DBTableRow::DBTableRow()
{

}
//____________________________________________________________________________
DBTableRow::DBTableRow(const DBTableFields * fields, TSQLRow * row) :
_fields(fields)
{
  for(unsigned int ifield = 0;
                 ifield < _fields->NFields(); ifield++)
                                        _row.push_back(row->GetField(ifield));
}
//____________________________________________________________________________
DBTableRow::DBTableRow(const DBTableRow * db_row)
{
  vector<string>::const_iterator row_iter;

  for(row_iter = db_row->_row.begin();
                  row_iter != db_row->_row.end(); ++row_iter)
                                            this->_row.push_back( *row_iter );
  this->_fields = db_row->_fields;
}
//____________________________________________________________________________  
DBTableRow::~DBTableRow()
{
  delete _fields;
}
//____________________________________________________________________________
string DBTableRow::Field(unsigned int field_pos) const
{
  return _row[field_pos];
}
//____________________________________________________________________________
string DBTableRow::Field(string field_name) const
{
  unsigned int field_pos = _fields->FieldPos(field_name);

  return Field(field_pos);
}
//____________________________________________________________________________
void DBTableRow::Print(ostream & stream) const
{
  vector<string>::const_iterator row_iter;

  for(row_iter = _row.begin(); row_iter != _row.end(); ++row_iter) {
       stream << *row_iter << endl;
  }  
}
//____________________________________________________________________________


