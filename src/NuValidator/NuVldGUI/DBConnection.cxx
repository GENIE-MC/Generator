//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBConnection

\brief    

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include "NuVldGUI/DBConnection.h"

using namespace genie::nuvld;

ClassImp(DBConnection)

//_____________________________________________________________________________
DBConnection::DBConnection()
{
  this->Init();
}
//_____________________________________________________________________________
DBConnection::~DBConnection()
{

}
//_____________________________________________________________________________
string DBConnection::URL(void) const
{
  ostringstream url;

  url << "mysql://" << _host << "/" << _database;

  return url.str();
}
//_____________________________________________________________________________
void DBConnection::Init(void)
{
  _host        = "";
  _database    = "";
  _user        = "";
  _password    = "";
  _sql_server  = 0;
}
//_____________________________________________________________________________

