//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Igor Kakorin <kakorin@jinr.ru>
 Joint Institute for Nuclear Research
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Physics/Common/NormXSec.h"

using namespace genie;
//____________________________________________________________________________
NormXSec::NormXSec() :
XSecAlgorithmI("genie::NormXSec")
{

}
//____________________________________________________________________________
NormXSec::NormXSec(string config) :
XSecAlgorithmI("genie::NormXSec", config)
{

}
//____________________________________________________________________________
NormXSec::~NormXSec()
{

}
//____________________________________________________________________________
double NormXSec::XSec(
   const Interaction * interaction, KinePhaseSpace_t kps) const
{
  const Target & tgt = interaction->InitState().Tgt();
  double A = tgt.A();
  return A*fNormScale*1e-11;
}
//_________________________________
double NormXSec::Integral(const Interaction * interaction) const
{
  const Target & tgt = interaction->InitState().Tgt();
  double A = tgt.A();
  return A*fNormScale*1e-11;
}
//____________________________________________________________________________
bool NormXSec::ValidProcess(const Interaction * interaction) const
{
  return true;
}
//____________________________________________________________________________
bool NormXSec::ValidKinematics(const Interaction * interaction) const
{
  return true;
}
//____________________________________________________________________________
void NormXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NormXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NormXSec::LoadConfig(void)
{
  GetParamDef( "NormScale", fNormScale, 1.0);
}
