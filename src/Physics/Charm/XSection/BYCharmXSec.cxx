#include <TMath.h>
#include <Math/Integrator.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/Charm/XSection/BYCharmXSec.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
/////#include "Numerical/IntegratorI.h"
#include "Physics/PartonDistributions/PDF.h"
#include "Physics/PartonDistributions/PDFModelI.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/DeepInelastic/XSection/QPMDISStrucFuncBase.h"

using namespace genie;
using namespace genie::constants;


//____________________________________________________________________________

double BYCharmXSec::XSec(const Interaction * interaction, KinePhaseSpace_t kps) const {
  if (!this->ValidProcess(interaction) || !this->ValidKinematics(interaction)){
    return 0;
  }
  const InitialState & init_state = interaction->InitState();
  const Target & tgt = init_state.Tgt();
  int  qpdg = tgt.HitQrkPdg();
  if (!pdg::IsCQuark(qpdg) &&  !pdg::IsAntiCQuark(qpdg)){
    return 0;
  }
  return this->fDISModel->XSec(interaction, kps);

}

void BYCharmXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BYCharmXSec::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
void BYCharmXSec::LoadConfig(void)
{
  fDISModel = 0;
  fDISModel =
      dynamic_cast<const XSecAlgorithmI *>(this->SubAlg("DISAlg"));
  assert(fDISModel);
}
