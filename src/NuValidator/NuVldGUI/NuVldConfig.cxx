//_____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldConfig

\brief    NuValidator GUI configuration options

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  September 17, 2005
*/
//_____________________________________________________________________________

#include <TSystem.h>

#include "Messenger/Messenger.h"
#include "NuVldGUI/NuVldConfig.h"

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

  string facade_enabled  = string( gSystem->Getenv("GOPT_ENABLE_NEUGEN") );

  LOG("NuVld", pDEBUG) << "GOPT_ENABLE_NEUGEN: " << facade_enabled;

  bool activate = ( facade_enabled.find("YES") != string::npos );

  if (!activate) {
     LOG("NuVld", pNOTICE) << "NeuGEN will be disabled in the NuVld session";
  }

  this->SetUseNeuGEN(activate);
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
