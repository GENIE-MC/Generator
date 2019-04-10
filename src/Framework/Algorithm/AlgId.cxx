//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 20, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include "Framework/Algorithm/AlgId.h"

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
AlgId::AlgId(const RgAlg & registry_item)
{
  this->Copy(registry_item);
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

  this->UpdateKey();
}
//____________________________________________________________________________
void AlgId::Copy(const RgAlg & registry_item)
{
  this->fName   = registry_item.name;
  this->fConfig = registry_item.config;

  this->UpdateKey();
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
