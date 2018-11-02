//_________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.
*/
//_________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/HadronTensors/HadronTensorPool.h"
#include "Physics/HadronTensors/ValenciaHadronTensorI.h"
#include "Physics/Multinucleon/XSection/SuSAv2MECPXSec.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

using namespace genie;

//_________________________________________________________________________
SuSAv2MECPXSec::SuSAv2MECPXSec() : XSecAlgorithmI("genie::SuSAv2MECPXSec")
{
}
//_________________________________________________________________________
SuSAv2MECPXSec::SuSAv2MECPXSec(string config)
  : XSecAlgorithmI("genie::SuSAv2MECPXSec", config)
{
}
//_________________________________________________________________________
SuSAv2MECPXSec::~SuSAv2MECPXSec()
{
}
//_________________________________________________________________________
double SuSAv2MECPXSec::XSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{
  // Get the hadron tensor pool
  HadronTensorPool& htp = HadronTensorPool::Instance();

  // Get the hadron tensor for the selected nuclide. Check the probe PDG code
  // to know whether to use the tensor for CC neutrino scattering or for
  // electron scattering
  int target_pdg = interaction->InitState().Tgt().Pdg();
  int probe_pdg = interaction->InitState().ProbePdg();

  HadronTensorType_t tensor_type = kHT_Undefined;
  if ( pdg::IsNeutrino(probe_pdg) || pdg::IsAntiNeutrino(probe_pdg) ) {
    tensor_type = kHT_MEC_FullAll;
  }
  else {
    // If the probe is not a neutrino, assume that it's an electron
    tensor_type = kHT_MEC_EM;
  }

  // The SuSAv2-MEC hadron tensors are defined using the same conventions
  // as the Valencia MEC model, so we can use the same sort of tensor
  // object to describe them.
  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>( htp.GetTensor(target_pdg,
    tensor_type, fHadronTensorTableName) );

  /// \todo Add scaling of the hadron tensor elements if an exact match
  /// for the requested target PDG code cannot be found in the hadron tensor
  /// table for this model.

  // If retrieving the tensor failed, complain and return zero
  if ( !tensor ) {
    LOG("SuSAv2MEC", pWARN) << "Failed to load a hadronic tensor for the"
      " nuclide " << target_pdg;
    return 0.;
  }

  // Check that the input kinematical point is within the range
  // in which hadron tensors are known (for chosen target)
  double Ev    = interaction->InitState().ProbeE(kRfLab);
  double Tl    = interaction->Kine().GetKV(kKVTl);
  double costl = interaction->Kine().GetKV(kKVctl);
  double ml    = interaction->FSPrimLepton()->Mass();
  double Q0    = 0.;
  double Q3    = 0.;

  genie::utils::mec::Getq0q3FromTlCostl(Tl, costl, Ev, ml, Q0, Q3);

  double Q0min = tensor->q0Min();
  double Q0max = tensor->q0Max();
  double Q3min = tensor->qMagMin();
  double Q3max = tensor->qMagMax();
  if (Q0 < Q0min || Q0 > Q0max || Q3 < Q3min || Q3 > Q3max) {
    return 0.0;
  }

  // Get the Q-value needed to calculate the cross sections using the
  // hadron tensor.
  /// \todo Shouldn't we get this from the nuclear model?
  int nu_pdg = interaction->InitState().ProbePdg();
  double Q_value = genie::utils::mec::Qvalue(target_pdg, nu_pdg);

  // Compute the cross section using the hadron tensor
  double xsec = tensor->dSigma_dT_dCosTheta(interaction, Q_value);

  // Apply the squared CKM matrix element Vud as a correction factor
  // (not included in the tabulated tensor element values)
  xsec *= fVud2;

  // Apply given overall scaling factor
  xsec *= fXSecScale;

  if ( kps != kPSTlctl ) {
    LOG("SuSAv2MEC", pWARN)
      << "Doesn't support transformation from "
      << KinePhaseSpace::AsString(kPSTlctl) << " to "
      << KinePhaseSpace::AsString(kps);
    xsec = 0.;
  }

  return xsec;
}
//_________________________________________________________________________
double SuSAv2MECPXSec::Integral(const Interaction* interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//_________________________________________________________________________
bool SuSAv2MECPXSec::ValidProcess(const Interaction* interaction) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const ProcessInfo& proc_info = interaction->ProcInfo();
  if ( !proc_info.IsMEC() ) {
    return false;
  }

  /// \todo Check whether CC, NC, EM? No tensor files for NC yet.

  return true;
}
//_________________________________________________________________________
void SuSAv2MECPXSec::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SuSAv2MECPXSec::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//_________________________________________________________________________
void SuSAv2MECPXSec::LoadConfig(void)
{
  // CKM matrix element connecting the up and down quarks
  // (not included in the tabulated SuSAv2-MEC hadron tensors)
  double Vud;
  GetParam("CKM-Vud", Vud);
  fVud2 = std::pow(Vud, 2);

  // Cross section scaling factor
  GetParam("MEC-XSecScale", fXSecScale) ;

  // Name of the hadron tensor table to use for this model
  GetParamDef("HadronTensorTableName", fHadronTensorTableName,
    std::string("SuSAv2_MEC"));

  fXSecIntegrator = dynamic_cast<const XSecIntegratorI*> (
    this->SubAlg("NumericalIntegrationAlg"));
  assert(fXSecIntegrator);
}
