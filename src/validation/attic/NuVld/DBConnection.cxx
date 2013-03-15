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

#include  "ValidationTools/NuVld/DBConnection.h"

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

