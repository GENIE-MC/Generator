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

#include <cassert>

#include "ValidationTools/NuVld/DBTableFields.h"

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
