//_____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldConfig

\brief    NuValidator GUI configuration options

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

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
  // if you can not find the NeuGEN library and NuValidator's NeuGEN,
  // deactivate NeuGEN

  string genie_path  = string( gSystem->Getenv("GENIE")       );
  string neugen_path = string( gSystem->Getenv("NEUGEN3PATH") );

  string nuvldn_lib = genie_path  + string("/lib/libGNuVldNeugen.so");
  string neugen_lib = neugen_path + string("/lib/libneugen3.a");

  LOG("NuVld", pDEBUG) << "NuVld/NeuGEN lib.: " << nuvldn_lib;
  LOG("NuVld", pDEBUG) << "NeuGEN lib.......: " << neugen_lib;

  bool nuvldn_lib_exists = ! (gSystem->AccessPathName( nuvldn_lib.c_str() ));
  bool neugen_lib_exists = ! (gSystem->AccessPathName( neugen_lib.c_str() ));

  bool activate = nuvldn_lib_exists && neugen_lib_exists;

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
