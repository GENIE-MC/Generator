#include "Tools/VMC/EvtLibPXSec.h"

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Physics/Common/VertexGenerator.h"
#include "Framework/GHEP/GHepParticle.h"

#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"

using namespace genie;
using namespace genie::vmc;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
EvtLibPXSec::EvtLibPXSec() :
XSecAlgorithmI("genie::vmc::EvtLibPXSec")
{

}
//____________________________________________________________________________
EvtLibPXSec::EvtLibPXSec(string config) :
XSecAlgorithmI("genie::vmc::EvtLibPXSec", config)
{

}
//____________________________________________________________________________
EvtLibPXSec::~EvtLibPXSec()
{

}
//____________________________________________________________________________
double EvtLibPXSec::XSec(const Interaction* in, KinePhaseSpace_t /*kps*/) const
{
  const InitialState& init_state = in->InitState();
  const double E  = init_state.ProbeE(kRfHitNucRest);
  const Target& target = init_state.Tgt();
  return target.N() * E; // TODO units
}

//____________________________________________________________________________
double EvtLibPXSec::Integral(const Interaction* in) const
{
  const InitialState& init_state = in->InitState();
  const double E  = init_state.ProbeE(kRfHitNucRest);
  const Target& target = init_state.Tgt();
  return target.N() * E*E/2; // TODO units
}
//____________________________________________________________________________
bool EvtLibPXSec::ValidProcess(const Interaction* /*in*/) const
{
  return true;
}
//____________________________________________________________________________
void EvtLibPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
}
//____________________________________________________________________________
void EvtLibPXSec::Configure(string config)
{
  Algorithm::Configure(config);
}
