//_________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.
*/
//_________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Multinucleon/XSection/NievesSimoVacasMECPXSec2016.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Physics/Multinucleon/XSection/HadronTensorI.h"
#include "Physics/Multinucleon/XSection/HadronTensorPool.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//_________________________________________________________________________
NievesSimoVacasMECPXSec2016::NievesSimoVacasMECPXSec2016() :
XSecAlgorithmI("genie::NievesSimoVacasMECPXSec2016")
{

}
//_________________________________________________________________________
NievesSimoVacasMECPXSec2016::NievesSimoVacasMECPXSec2016(string config) :
XSecAlgorithmI("genie::NievesSimoVacasMECPXSec2016", config)
{

}
//_________________________________________________________________________
NievesSimoVacasMECPXSec2016::~NievesSimoVacasMECPXSec2016()
{

}
//_________________________________________________________________________
double NievesSimoVacasMECPXSec2016::XSec(
  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  // This function returns d2sigma/(dTmu dcos_mu) in GeV^(-3)

  // Get the hadron tensor pool
  HadronTensorPool& htp = HadronTensorPool::Instance();

  int target_pdg = interaction->InitState().Tgt().Pdg();

  int A_request = pdg::IonPdgCodeToA(target_pdg);
  int Z_request = pdg::IonPdgCodeToZ(target_pdg);

  // To generate cross-sections for nuclei other than those with hadron
  // tensors we need to pull both the full cross-section and
  // the pn initial state fraction.
  // Non-isoscalar nuclei are beyond the original published Valencia model
  // and scale with A according to the number of pp, pn, or nn pairs
  // the probe is expected to find.
  // There is some by-hand optimization here, skipping the delta part when
  // only the total cross-section is requested.
  // Possible future models without a Delta had tensor would also use that
  // flag to call this without computing the Delta part.

  // Try to look up a hadron tensor in the pool that is an exact match for
  // the target nucleus. If an exact match cannot be found, decide upon a
  // suitable substitute based on the mass number A and proton number Z.

  int tensor_pdg = target_pdg;

  /// \todo Replace these hard-coded replacements with an equivalent XML
  /// configuration
  if ( ! htp.GetTensor(tensor_pdg, HadronTensorType::kHT_MEC_FullAll) ) {

    if ( A_request == 4 && Z_request == 2 ) {
      tensor_pdg = kPdgTgtC12;
      // This is for helium 4, but use carbon tensor
      // the use of nuclear density parameterization is suspicious
      // but some users (MINERvA) need something not nothing.
      // The pn will be exactly 1/3, but pp and nn will be ~1/4
      // Because the combinatorics are different.
      // Could do lithium beryllium boron which you don't need
    }
    else if (A_request < 9) {
      // refuse to do D, T, He3, Li, and some Be, B
      // actually it would work technically, maybe except D, T
      MAXLOG("NievesSimoVacasMEC", pWARN, 10)
          << "Asked to scale to deuterium through boron "
          << target_pdg << " nope, lets not do that.";
      return 0;
    }
    else if (A_request >= 9 && A_request < 15) {
      tensor_pdg = kPdgTgtC12;
      //}
      // could explicitly put in nitrogen for air
      //else if ( A_request >= 14 && A < 15) { // AND CHANGE <=14 to <14.
      //  tensor_pdg = kPdgTgtN14;
    }
    else if (A_request >= 15 && A_request < 22) {
      tensor_pdg = kPdgTgtO16;
    }
    else if (A_request >= 22 && A_request < 33) {
      // of special interest, this gets Al27 and Si28
      tensor_pdg = 1000140280;
    }
    else if (A_request >= 33 && A_request < 50) {
      // of special interest, this gets Ar40 and Ti48
      tensor_pdg = kPdgTgtCa40;
    }
    else if (A_request >= 50 && A_request < 90) {
      // pseudoFe56, also covers many other ferrometals and Ge
      tensor_pdg = 1000280560;
    }
    else if (A_request >= 90 && A_request < 160) {
      // use Ba112 = PseudoCd.  Row5 of Periodic table useless. Ag, Xe?
      tensor_pdg = 1000561120;
    }
    else if (A_request >= 160) {
      // use Rf208 = pseudoPb
      tensor_pdg = 1001042080;
    }
    else {
      MAXLOG("NievesSimoVacasMEC", pWARN, 10)
          << "Asked to scale to a nucleus "
          << target_pdg << " which we don't know yet.";
      return 0;
    }
  }

  // Check that the input kinematical point is within the range
  // in which hadron tensors are known (for chosen target)
  double Ev    = interaction->InitState().ProbeE(kRfLab);
  double Tl    = interaction->Kine().GetKV(kKVTl);
  double costl = interaction->Kine().GetKV(kKVctl);
  double ml    = interaction->FSPrimLepton()->Mass();
  double Q0    = 0;
  double Q3    = 0;

  genie::utils::mec::Getq0q3FromTlCostl(Tl, costl, Ev, ml, Q0, Q3);

  const HadronTensorI* tensor = htp.GetTensor(tensor_pdg,
    HadronTensorType::kHT_MEC_FullAll);

  // If retrieving the tensor failed, complain and return zero
  if ( !tensor ) {
    LOG("NievesSimoVacasMEC", pWARN) << "Failed to load a"
      " hadronic tensor for the nuclide " << tensor_pdg;
    return 0.;
  }

  // Assume for now that the range of validity for the "FullAll" hadron
  // tensor is the same as for the partial hadron tensors
  /// \todo Revisit this assumption, and perhaps implement something
  /// more robust
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

  // By default, we will compute the full cross-section. If a resonance is
  // set, we will calculate the part of the cross-section with an internal
  // Delta line without a final state pion (usually called PPD for pioness
  // Delta decay). If a {p,n} hit dinucleon was set we will calculate the
  // cross-section for that component only (either full or PDD cross-section)
  bool delta = interaction->ExclTag().KnownResonance();
  bool pn    = (interaction->InitState().Tgt().HitNucPdg() == kPdgClusterNP);

  double xsec_all = 0.;
  double xsec_pn  = 0.;
  if ( delta ) {

    const HadronTensorI* tensor_delta_all = htp.GetTensor(tensor_pdg,
      HadronTensorType::kHT_MEC_DeltaAll);

    if ( !tensor_delta_all ) {
      LOG("NievesSimoVacasMEC", pWARN) << "Failed to load a \"DeltaAll\""
        << " hadronic tensor for nuclide " << tensor_pdg;
      return 0.;
    }

    const HadronTensorI* tensor_delta_pn = htp.GetTensor(tensor_pdg,
      HadronTensorType::kHT_MEC_Deltapn);

    if ( !tensor_delta_pn ) {
      LOG("NievesSimoVacasMEC", pWARN) << "Failed to load a \"Deltapn\""
        << " hadronic tensor for nuclide " << tensor_pdg;
      return 0.;
    }

    xsec_all = tensor_delta_all->dSigma_dT_dCosTheta(interaction, Q_value);
    xsec_pn = tensor_delta_pn->dSigma_dT_dCosTheta(interaction, Q_value);

  }
  else {

    const HadronTensorI* tensor_full_all = htp.GetTensor(tensor_pdg,
      HadronTensorType::kHT_MEC_FullAll);

    if ( !tensor_full_all ) {
      LOG("NievesSimoVacasMEC", pWARN) << "Failed to load a \"FullAll\""
        << " hadronic tensor for nuclide " << tensor_pdg;
      return 0.;
    }

    const HadronTensorI* tensor_full_pn = htp.GetTensor(tensor_pdg,
      HadronTensorType::kHT_MEC_Fullpn);

    if ( !tensor_full_pn ) {
      LOG("NievesSimoVacasMEC", pWARN) << "Failed to load a \"Fullpn\""
        << " hadronic tensor for nuclide " << tensor_pdg;
      return 0.;
    }

    xsec_all = tensor_full_all->dSigma_dT_dCosTheta(interaction, Q_value);
    xsec_pn = tensor_full_pn->dSigma_dT_dCosTheta(interaction, Q_value);
  }

  // We need to scale the cross section appropriately if
  // we are using a hadronic tensor for a nuclide that is different
  // from the actual target
  bool need_to_scale = (target_pdg != tensor_pdg);

  // would need to trap and treat He3, T, D special here.
  if ( need_to_scale ) {

    double PP = Z_request;
    double NN = A_request - PP;
    double P  = pdg::IonPdgCodeToZ(tensor_pdg);
    double N  = pdg::IonPdgCodeToA(tensor_pdg) - P;

    double scale_pn = TMath::Sqrt( (PP*NN)/(P*N) );
    double scale_pp = TMath::Sqrt( (PP * (PP - 1.)) / (P * (P - 1.)) );
    double scale_nn = TMath::Sqrt( (NN * (NN - 1.)) / (N * (N - 1.)) );

    LOG("NievesSimoVacasMEC", pDEBUG)
        << "Scale pn pp nn for (" << target_pdg << ", " << tensor_pdg << ")"
        << " : " << scale_pn << " " << scale_pp << " " << scale_nn;

    // This is an approximation in at least three senses:
    // 1. We are scaling from an isoscalar nucleus using p and n counting
    // 2. We are not using the right qvalue in the had tensor
    // 3. We are not scaling the Delta faster than the non-Delta.
    // The guess is that these are good approximations.
    // A test we could document is to scale from O16 to N14 or C12 using this
    // algorithm and see how many percent deviation we see from the full
    // calculation.
    double temp_all = xsec_all;
    double temp_pn  = xsec_pn * scale_pn;
    if (nu_pdg > 0) {
      // matter neutrinos
      temp_all = xsec_pn * scale_pn + (xsec_all - xsec_pn) * scale_nn;
    }
    else {
      // antineutrinos
      temp_all = xsec_pn * scale_pn + (xsec_all - xsec_pn) * scale_pp;
    }

    xsec_all = temp_all;
    xsec_pn  = temp_pn;

  }

  // Choose the right kind of cross section ("all" or "pn") to return
  // based on whether a {p, n} dinucleon was hit
  double xsec = (pn) ? xsec_pn : xsec_all;

  // Apply given scaling factor
  xsec *= fXSecScale;

  if ( kps != kPSTlctl ) {
    LOG("NievesSimoVacasMEC", pWARN)
      << "Doesn't support transformation from "
      << KinePhaseSpace::AsString(kPSTlctl) << " to "
      << KinePhaseSpace::AsString(kps);
    xsec = 0;
  }

  return xsec;
}
//_________________________________________________________________________
double NievesSimoVacasMECPXSec2016::Integral(
        const Interaction * interaction) const
{
    double xsec = fXSecIntegrator->Integrate(this,interaction);
    return xsec;
}
//_________________________________________________________________________
bool NievesSimoVacasMECPXSec2016::ValidProcess(
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
void NievesSimoVacasMECPXSec2016::Configure(const Registry & config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void NievesSimoVacasMECPXSec2016::Configure(string config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//_________________________________________________________________________
void NievesSimoVacasMECPXSec2016::LoadConfig(void)
{
	// Cross section scaling factor
	GetParam( "MEC-CC-XSecScale", fXSecScale ) ;

	fXSecIntegrator =
        dynamic_cast<const XSecIntegratorI *> (
                this->SubAlg("NumericalIntegrationAlg"));
    assert(fXSecIntegrator);

}
//_________________________________________________________________________
