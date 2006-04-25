//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 20, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include "Algorithm/AlgId.h"

using std::ostringstream;

using namespace genie;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const AlgId & algid)
  {
    algid.Print(stream);
    return stream;
  }
}
//____________________________________________________________________________
AlgId::AlgId()
{
  this->Init();
}
//____________________________________________________________________________
AlgId::AlgId(string name, string config)
{
  this->SetId(name,config);
}
//____________________________________________________________________________
AlgId::AlgId(const AlgId & id)
{
  this->Copy(id);
}
//____________________________________________________________________________
AlgId::~AlgId()
{

}
//____________________________________________________________________________
void AlgId::SetName(string name)
{
  this->SetId(name, this->Config());
}
//____________________________________________________________________________
void AlgId::SetConfig(string config)
{
  this->SetId(this->Name(), config);
}
//____________________________________________________________________________
void AlgId::SetId(string name, string config)
{
  this->fName   = name;
  this->fConfig = config;

  this->UpdateKey();
}
//____________________________________________________________________________
void AlgId::Copy(const AlgId & id)
{
  this->fName   = id.Name();
  this->fConfig = id.Config();
  this->fKey    = id.Key();
}
//____________________________________________________________________________
void AlgId::Print(ostream & stream) const
{
  stream << this->Key();
}
//____________________________________________________________________________
void AlgId::UpdateKey(void)
{
  ostringstream key;

  key << this->Name();
  if(this->Config().size() > 0) key << "/" << this->Config();

  fKey = key.str();
}
//____________________________________________________________________________
void AlgId::Init(void)
{
  this->fName   = "";
  this->fConfig = "";
  this->fKey    = "";
}
//____________________________________________________________________________
