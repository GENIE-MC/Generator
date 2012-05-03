//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBQueryString

\brief    A collectable string class with encoded information needed to form
          a NuVld database SQL query

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 26, 2005
*/
//_____________________________________________________________________________

#ifndef _DB_QUERY_STRING_H_
#define _DB_QUERY_STRING_H_

#include <iostream>
#include <string>

#include <TObjString.h>

#include  "ValidationTools/NuVld/DBTableType.h"

using std::ostream;
using std::string;

namespace genie {
namespace nuvld {

class DBQueryString : public TObjString
{
public:

   DBQueryString(void);
   DBQueryString(const DBQueryString & query_string);
   DBQueryString(const char * query_string);
   DBQueryString(string query_string);
   virtual ~DBQueryString();

   bool          IsNull    (void) const;
   string        AsString  (void) const;  
   string        KeyList   (void) const;
   string        CutList   (void) const;   
   DBTableType_t TableType (void) const;

   friend ostream & operator << (ostream & stream, const DBQueryString & query_string);
             
private:

   string ExtractFieldValue (string field_name) const;

ClassDef(DBQueryString, 1)  
};

} // nuvld namespace
} // genie namespace

#endif  // _DB_QUERY_STRING_H_
