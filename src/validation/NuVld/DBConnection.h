//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBConnection

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _DBASE_CONNECTION_H_
#define _DBASE_CONNECTION_H_

#include <string>
#include <sstream>

#include <TSQLServer.h>
#include <TObject.h>

using std::string;
using std::ostringstream;

namespace genie {
namespace nuvld {

class DBConnection : public TObject {

public:

   friend class DBConnectionDialog;
   
   DBConnection();
   virtual ~DBConnection();

   string       URL        (void) const;
   string       Host       (void) const { return _host;       }
   string       DataBase   (void) const { return _database;   }
   string       User       (void) const { return _user;       }
   string       Password   (void) const { return _password;   }
   TSQLServer * SqlServer  (void) const { return _sql_server; }

   void SetSqlServer (TSQLServer * db)  { _sql_server = db;  }

private:

   void Init (void);

   string       _host;
   string       _database;
   string       _user;
   string       _password;
   TSQLServer * _sql_server;

   ClassDef(DBConnection, 0);
};

} // nuvld namespace
} // genie namespace

#endif

