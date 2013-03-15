//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "ValidationTools/NuVld/DBTableRow.h"

using std::endl;
using namespace genie::nuvld;

ClassImp(DBTableRow)

//____________________________________________________________________________
DBTableRow::DBTableRow()
{

}
//____________________________________________________________________________
DBTableRow::DBTableRow(const DBTableFields * fields, TSQLRow * row) :
fFields(fields)
{
  for(unsigned int ifield = 0;
                 ifield < fFields->NFields(); ifield++)
                                        fRow.push_back(row->GetField(ifield));
}
//____________________________________________________________________________
DBTableRow::DBTableRow(const DBTableRow * db_row)
{
  vector<string>::const_iterator row_iter;

  for(row_iter = db_row->fRow.begin();
                  row_iter != db_row->fRow.end(); ++row_iter)
                                            this->fRow.push_back( *row_iter );
  this->fFields = db_row->fFields;
}
//____________________________________________________________________________
DBTableRow::~DBTableRow()
{
  delete fFields;
}
//____________________________________________________________________________
string DBTableRow::Field(unsigned int field_pos) const
{
  return fRow[field_pos];
}
//____________________________________________________________________________
string DBTableRow::Field(string field_name) const
{
  unsigned int field_pos = fFields->FieldPos(field_name);

  return Field(field_pos);
}
//____________________________________________________________________________
void DBTableRow::Print(ostream & stream) const
{
  vector<string>::const_iterator row_iter;

  for(row_iter = fRow.begin(); row_iter != fRow.end(); ++row_iter) {
       stream << *row_iter << endl;
  }
}
//____________________________________________________________________________


