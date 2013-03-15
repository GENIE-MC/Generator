//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Sep 17, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include <TSystem.h>

#include "Messenger/Messenger.h"
#include  "ValidationTools/NuVld/NuVldConfig.h"

using namespace genie::nuvld;

//_____________________________________________________________________________
NuVldConfig::NuVldConfig()
{
  this->Init();
}
//_____________________________________________________________________________
NuVldConfig::NuVldConfig(const NuVldConfig & config)
{
  this->Copy(config);
}
//_____________________________________________________________________________
NuVldConfig::~NuVldConfig()
{

}
//_____________________________________________________________________________
void NuVldConfig::AutoDetect(void)
{
// if you can not find the GENIE's NeuGEN facade was not enabled (and therefore
// a dummy NeuGEN library was built in the program) deactivate NeuGEN

/*
  string facade_enabled  = string( gSystem->Getenv("GOPT_ENABLE_NEUGEN") );

  LOG("NuVld", pDEBUG) << "GOPT_ENABLE_NEUGEN: " << facade_enabled;

  bool activate = ( facade_enabled.find("YES") != string::npos );

  if (!activate) {
     LOG("NuVld", pNOTICE) << "NeuGEN will be disabled in the NuVld session";
  }

  this->SetUseNeuGEN(activate);
*/

  this->SetUseNeuGEN(false);
}
//_____________________________________________________________________________
void NuVldConfig::Init(void)
{
  this->SetUseNeuGEN        (true);
  this->SetUseCompactLayout (false);
}
//_____________________________________________________________________________
void NuVldConfig::Copy(const NuVldConfig & config)
{
  this->SetUseNeuGEN        (config.UseNeuGEN()        );
  this->SetUseCompactLayout (config.UseCompactLayout() );
}
//_____________________________________________________________________________
