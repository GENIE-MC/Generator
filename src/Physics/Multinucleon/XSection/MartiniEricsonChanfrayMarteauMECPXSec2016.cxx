//_________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.
*/
//_________________________________________________________________________

#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Multinucleon/XSection/MartiniEricsonChanfrayMarteauMECPXSec2016.h"
#include "Physics/Multinucleon/XSection/MECHadronTensor.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//_________________________________________________________________________
MartiniEricsonChanfrayMarteauMECPXSec2016::
  MartiniEricsonChanfrayMarteauMECPXSec2016() :
XSecAlgorithmI("genie::MartiniEricsonChanfrayMarteauMECPXSec2016")
{

}
//_________________________________________________________________________
MartiniEricsonChanfrayMarteauMECPXSec2016::
  MartiniEricsonChanfrayMarteauMECPXSec2016(string config) :
XSecAlgorithmI("genie::NievesSimoVacasMECPXSec2016",config)
{

}
//_________________________________________________________________________
MartiniEricsonChanfrayMarteauMECPXSec2016::
  ~MartiniEricsonChanfrayMarteauMECPXSec2016()
{

}
//_________________________________________________________________________
double MartiniEricsonChanfrayMarteauMECPXSec2016::XSec(
  const Interaction * /*interaction*/, KinePhaseSpace_t /*kps*/) const
{
// The basic quantity calculated is d2sigma/(dT_l dcostheta_l)

  return 0;
}
//_________________________________________________________________________
double MartiniEricsonChanfrayMarteauMECPXSec2016::Integral(
  const Interaction * interaction) const 
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//_________________________________________________________________________
bool MartiniEricsonChanfrayMarteauMECPXSec2016::ValidProcess(
  const Interaction * interaction) const 
{
  if (interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info = interaction->ProcInfo();
  if (!proc_info.IsMEC()) {
      return false;
  }
  return true;
}
//_________________________________________________________________________
void MartiniEricsonChanfrayMarteauMECPXSec2016::Configure(
  const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MartiniEricsonChanfrayMarteauMECPXSec2016::Configure(
  string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//_________________________________________________________________________
void MartiniEricsonChanfrayMarteauMECPXSec2016::LoadConfig(void)
{
  fXSecIntegrator =
     dynamic_cast<const XSecIntegratorI *> (
           this->SubAlg("NumericalIntegrationAlg"));
  assert(fXSecIntegrator);
}
//_________________________________________________________________________
