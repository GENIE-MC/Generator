//_________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
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
#include "Physics/HadronTensors/MartiniMECHadronTensorModel.h"
#include "Physics/HadronTensors/LabFrameHadronTensorI.h"
#include "Physics/Multinucleon/XSection/MartiniEricsonChanfrayMarteauMECPXSec2024.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

using namespace genie;


//_________________________________________________________________________
MartiniEricsonChanfrayMarteauMECPXSec2024::MartiniEricsonChanfrayMarteauMECPXSec2024() : XSecAlgorithmI("genie::MartiniEricsonChanfrayMarteauMECPXSec2024")
{
}
//_________________________________________________________________________
MartiniEricsonChanfrayMarteauMECPXSec2024::MartiniEricsonChanfrayMarteauMECPXSec2024(string config)
  : XSecAlgorithmI("genie::MartiniEricsonChanfrayMarteauMECPXSec2024", config)
{
}
//_________________________________________________________________________
MartiniEricsonChanfrayMarteauMECPXSec2024::~MartiniEricsonChanfrayMarteauMECPXSec2024()
{
}
//_________________________________________________________________________
double MartiniEricsonChanfrayMarteauMECPXSec2024::XSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{
  // Don't try to do the calculation if we've been handed an interaction that
  // doesn't make sense
  if ( !this->ValidProcess(interaction) ) return 0.;

  // Get the hadron tensor for the selected nuclide. Check the probe PDG code
  // to know whether to use the tensor for CC neutrino scattering or for
  // electron scattering
  int target_pdg = interaction->InitState().Tgt().Pdg();
  int probe_pdg = interaction->InitState().ProbePdg();
  int A_request = pdg::IonPdgCodeToA(target_pdg);
  int Z_request = pdg::IonPdgCodeToZ(target_pdg);
  bool need_to_scale = false;

  HadronTensorType_t tensor_type = kHT_Undefined;
  if ( pdg::IsNeutrino(probe_pdg) || pdg::IsAntiNeutrino(probe_pdg) ) {
    tensor_type = kHT_MEC_FullAll;
    //pn_tensor_type = kHT_MEC_Fullpn;
    //tensor_type = kHT_MEC_FullAll_wImag;
    //pn_tensor_type = kHT_MEC_FullAll_wImag;
  }
  else {
    // If the probe is not a neutrino, assume that it's an electron
    // For the moment all electron interactions are pp final state
    tensor_type = kHT_MEC_EM;
    //pn_tensor_type = kHT_MEC_EM;
  }

  int tensor_pdg = target_pdg;
  if ( ! fHadronTensorModel->GetTensor(tensor_pdg, genie::kHT_MEC_FullAll) )
  {

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
      MAXLOG("MartiniEricsonChanfrayMarteauMEC", pWARN, 10)
          << "Asked to scale to deuterium through boron "
          << target_pdg << " nope, lets not do that.";
      return 0;
    }
    else if (A_request >= 9 && A_request < 15) {
      tensor_pdg = kPdgTgtC12;
    }
    else if (A_request >= 15 && A_request < 22) {
      tensor_pdg = kPdgTgtO16;
    }
    else if (A_request >= 22 && A_request < 33) {
      // of special interest, this gets Al27 and Si28
      tensor_pdg = kPdgTgtO16;
    }
    else if (A_request >= 33 && A_request < 50) {
      // of special interest, this gets Ar40 and Ti48
      tensor_pdg = kPdgTgtCa40;
    }
    else if (A_request >= 50 && A_request < 90) {
      // pseudoFe56, also covers many other ferrometals and Ge
      tensor_pdg = kPdgTgtCa40;
    }
    else if (A_request >= 90 && A_request < 160) {
      // use Ba112 = PseudoCd.  Row5 of Periodic table useless. Ag, Xe?
      tensor_pdg = kPdgTgtCa40;
    }
    else if (A_request >= 160) {
      // use Rf208 = pseudoPb
      tensor_pdg = kPdgTgtCa40;
    }
    else {
      MAXLOG("MartiniEricsonChanfrayMarteauMEC", pWARN, 10)
          << "Asked to scale to a nucleus "
          << target_pdg << " which we don't know yet.";
      return 0;
    }
  }

  // The MartiniEricsonChanfrayMarteau-MEC hadron tensors are defined using the same conventions
  // as the Valencia MEC model, so we can use the same sort of tensor
  // object to describe them.
  const LabFrameHadronTensorI* tensor
    = dynamic_cast<const LabFrameHadronTensorI*>( fHadronTensorModel->GetTensor(tensor_pdg,
    tensor_type) );

  // If retrieving the tensor failed, complain and return zero
  if ( !tensor ) {
    LOG("MartiniEricsonChanfrayMarteauMEC", pWARN) << "Failed to load a hadronic tensor for the"
      " nuclide " << tensor_pdg;
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

  // The Q-Value essentially corrects q0 to account for nuclear
  // binding energy in the Valencia model but this effect is already
  // in Guille's tensors so its set it to 0.
  // However, additional corrections may be necessary:
  double Delta_Q_value = Qvalue( * interaction ) ;

  genie::utils::mec::Getq0q3FromTlCostl(Tl, costl, Ev, ml, Q0, Q3);

  double Q0min = tensor->q0Min();
  double Q0max = tensor->q0Max();
  double Q3min = tensor->qMagMin();
  double Q3max = tensor->qMagMax();
  if (Q0-Delta_Q_value < Q0min || Q0-Delta_Q_value > Q0max || Q3 < Q3min || Q3 > Q3max) {
    return 0.0;
  }

  // *** Enforce the global Q^2 cut (important for EM scattering) ***
  // Choose the appropriate minimum Q^2 value based on the interaction
  // mode (this is important for EM interactions since the differential
  // cross section blows up as Q^2 --> 0)
  double Q2min = genie::controls::kMinQ2Limit; // CC/NC limit
  if ( interaction->ProcInfo().IsEM() ) Q2min = genie::utils::kinematics
    ::electromagnetic::kMinQ2Limit; // EM limit

  // Neglect shift due to binding energy. The cut is on the actual
  // value of Q^2, not the effective one to use in the tensor contraction.
  double Q2 = Q3*Q3 - Q0*Q0;
  if ( Q2 < Q2min ) return 0.;

  // By default, we will compute the full cross-section. If a {p,n} hit
  // dinucleon was set we will calculate the cross-section for that
  // component only

  bool pn = (interaction->InitState().Tgt().HitNucPdg() == kPdgClusterNP);

  // Compute the cross section using the hadron tensor
  double xsec = tensor->dSigma_dT_dCosTheta_rosenbluth(interaction, Delta_Q_value);

  // This scaling should be okay-ish for the total xsec, but it misses
  // the energy shift. To get this we should really just build releveant
  // hadron tensors but there may be some ways to approximate it.
  // For more details see Guille's thesis: https://idus.us.es/xmlui/handle/11441/74826
  if ( need_to_scale ) {
    FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
    const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
    double KF_tgt = kft->FindClosestKF(target_pdg, kPdgProton);
    double KF_ten = kft->FindClosestKF(tensor_pdg, kPdgProton);
    LOG("MartiniEricsonChanfrayMarteauMEC", pDEBUG) << "KF_tgt = " << KF_tgt;
    LOG("MartiniEricsonChanfrayMarteauMEC", pDEBUG) << "KF_ten = " << KF_ten;
    double A_ten  = pdg::IonPdgCodeToA(tensor_pdg);
    double scaleFact = (A_request/A_ten)*(KF_tgt/KF_ten)*(KF_tgt/KF_ten);
    xsec *= scaleFact;
  }

  // Apply given overall scaling factor
  xsec *= fXSecScale;

  // Scale given a scaling algorithm:
  if( fMECScaleAlg ) xsec *= fMECScaleAlg->GetScaling( * interaction ) ;

  if ( kps != kPSTlctl ) {
    LOG("MartiniEricsonChanfrayMarteauMEC", pWARN)
      << "Doesn't support transformation from "
      << KinePhaseSpace::AsString(kPSTlctl) << " to "
      << KinePhaseSpace::AsString(kps);
    xsec = 0.;
  }

  return xsec;
}
//_________________________________________________________________________
double MartiniEricsonChanfrayMarteauMECPXSec2024::PairRatio(const Interaction* interaction,
  const std::string& final_state_ratio) const
{

  // Currently we only have the relative pair contributions for C12.
  // We hope to add mode later, but for the moment assume the relative
  // contributions are A-independant.

  int probe_pdg = interaction->InitState().ProbePdg();

  HadronTensorType_t tensor_type = kHT_Undefined;
  HadronTensorType_t pn_tensor_type = kHT_Undefined;
  HadronTensorType_t pp_tensor_type = kHT_Undefined;

  if ( pdg::IsNeutrino(probe_pdg) || pdg::IsAntiNeutrino(probe_pdg) ) {
    tensor_type = kHT_MEC_FullAll;
    pn_tensor_type = kHT_MEC_Fullpn;
  }
  else {
    // If the probe is not a neutrino, assume that it's an electron
    // For the moment all electron interactions are pp final state
    tensor_type = kHT_MEC_EM;
    pn_tensor_type = kHT_MEC_EM_pn;
    pp_tensor_type = kHT_MEC_EM_pp;
  }

  // The MartiniEricsonChanfrayMarteau-MEC hadron tensors are defined using the same conventions
  // as the Valencia MEC model, so we can use the same sort of tensor
  // object to describe them.
  const LabFrameHadronTensorI* tensor
    = dynamic_cast<const LabFrameHadronTensorI*>( fHadronTensorModel->GetTensor(kPdgTgtC12,
    tensor_type) );

  const LabFrameHadronTensorI* tensor_pn
    = dynamic_cast<const LabFrameHadronTensorI*>( fHadronTensorModel->GetTensor(kPdgTgtC12,
    pn_tensor_type) );

  const LabFrameHadronTensorI* tensor_pp
    = dynamic_cast<const LabFrameHadronTensorI*>( fHadronTensorModel->GetTensor(kPdgTgtC12,
    pp_tensor_type) );

  // If retrieving the tensor failed, complain and return zero
  if ( !tensor ) {
    LOG("MartiniEricsonChanfrayMarteauMEC", pWARN) << "Failed to load a hadronic tensor for the"
      " nuclide " << kPdgTgtC12;
    return 0.;
  }

  if ( !tensor_pn ) {
    LOG("MartiniEricsonChanfrayMarteauMEC", pWARN) << "Failed to load pn hadronic tensor for the"
      " nuclide " << kPdgTgtC12;
    return 0.;
  }

  if ( !tensor_pp && interaction->ProcInfo().IsEM() ) {
    LOG("MartiniEricsonChanfrayMarteauMEC", pWARN) << "Failed to load pp hadronic tensor for the"
      " nuclide " << kPdgTgtC12;
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
    return 1.0;
  }

  // The Q-Value essentially corrects q0 to account for nuclear
  // binding energy in the Valencia model but this effect is already
  // in Guille's tensors so its set it to 0.
  // However, additional corrections may be necessary:
  double Delta_Q_value = Qvalue( * interaction ) ;

  // Compute the cross section using the hadron tensor
  double xsec_all = tensor->dSigma_dT_dCosTheta_rosenbluth(interaction, Delta_Q_value);

  double ratio;

  if (final_state_ratio == "pnFraction") { // pnFraction will be calculated by default
    double xsec_pn = tensor_pn->dSigma_dT_dCosTheta_rosenbluth(interaction, Delta_Q_value);

    //hadron tensor precision can sometimes lead to 0 xsec_pn but finite xsec
    //seems to cause issues downstream ...
    if(xsec_pn==0) xsec_pn = 0.00001*xsec_all;

    double pn_ratio = (1e10*xsec_pn)/(1e10*xsec_all);

    ratio = pn_ratio;

  } else if (final_state_ratio == "ppFraction") {
    double xsec_pp = tensor_pp->dSigma_dT_dCosTheta_rosenbluth(interaction, Delta_Q_value);

    if(xsec_pp==0) xsec_pp = 0.00001*xsec_all;

    double pp_ratio = (1e10*xsec_pp)/(1e10*xsec_all);

    ratio = pp_ratio;

  }

  return ratio;
}
//_________________________________________________________________________
double MartiniEricsonChanfrayMarteauMECPXSec2024::Integral(const Interaction* interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//_________________________________________________________________________
bool MartiniEricsonChanfrayMarteauMECPXSec2024::ValidProcess(const Interaction* interaction) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const ProcessInfo& proc_info = interaction->ProcInfo();
  if ( !proc_info.IsMEC() ) {
    return false;
  }

  int probe = interaction->InitState().ProbePdg();

  bool is_nu = pdg::IsNeutrino( probe );
  bool is_nub = pdg::IsAntiNeutrino( probe );
  bool is_chgl = pdg::IsChargedLepton( probe );

  bool prc_ok = ( proc_info.IsWeakCC() && (is_nu || is_nub) )
    || ( proc_info.IsEM() && is_chgl );

  if ( !prc_ok ) return false;

  return true;
}
//_________________________________________________________________________
double MartiniEricsonChanfrayMarteauMECPXSec2024::Qvalue(const Interaction & interaction ) const
{
  // Get the hadron tensor for the selected nuclide. Check the probe PDG code
  // to know whether to use the tensor for CC neutrino scattering or for
  // electron scattering
  int target_pdg = interaction.InitState().Tgt().Pdg();
  int probe_pdg = interaction.InitState().ProbePdg();
  int tensor_pdg = kPdgTgtC12;
  int A_request = pdg::IonPdgCodeToA(target_pdg);

  double Eb_tgt=0;
  double Eb_ten=0;

  /// \todo Add more hadron tensors so this scaling is not so terrible
  // At the moment all we have is Carbon so this is all just a place holder ...
  if ( A_request == 4 ) {
    Eb_tgt=fEbHe; Eb_ten=fEbC;
    // This is for helium 4, but use carbon tensor, may not be ideal ...
  }
  else if (A_request < 9) {
    Eb_tgt=fEbLi; Eb_ten=fEbC;
  }
  else if (A_request >= 9 && A_request < 15) {
    Eb_tgt=fEbC; Eb_ten=fEbC;
  }
  else if(A_request >= 15 && A_request < 22) {
    Eb_tgt=fEbO; Eb_ten=fEbO;
  }
  else if(A_request >= 22 && A_request < 40) {
    Eb_tgt=fEbMg; Eb_ten=fEbO;
  }
  else if(A_request >= 40 && A_request < 56) {
    Eb_tgt=fEbCa; Eb_ten=fEbCa;
  }
  else if(A_request >= 56 && A_request < 119) {
    Eb_tgt=fEbFe; Eb_ten=fEbCa;
  }
  else if(A_request >= 119 && A_request < 206) {
    Eb_tgt=fEbSn; Eb_ten=fEbCa;
  }
  else if(A_request >= 206) {
    Eb_tgt=fEbPb; Eb_ten=fEbCa;
  }

  // SD: The Q-Value essentially corrects q0 to account for nuclear
  // binding energy in the Valencia model but this effect is already
  // in Guille's tensors so I'll set it to 0.
  // However, if I want to scale I need to account for the altered
  // binding energy. To first order I can use the Delta_Q_value for this.
  // But this is 2p2h - so binding energy counts twice - use 2*1p1h
  // value (although what should be done here is still not clear).

  double Delta_Q_value = 2*(Eb_tgt-Eb_ten);

  // Apply Qvalue relative shift if needed:
  if( fQvalueShifter ) {
    // We have the option to add an additional shift on top of the binding energy correction
    // The QvalueShifter, is a relative shift to the Q_value.
    // The Q_value was already taken into account in the hadron tensor. Here we recalculate it
    // to get the right absolute shift.
    double tensor_Q_value = genie::utils::mec::Qvalue(tensor_pdg,probe_pdg);
    double total_Q_value = tensor_Q_value + Delta_Q_value ;
    double Q_value_shift = total_Q_value * fQvalueShifter -> Shift( interaction.InitState().Tgt() ) ;
    Delta_Q_value += Q_value_shift ;
  }

  // We apply an extra Q-value shift here to account for differences between
  // the 12C EM MEC tensors currently in use (which have a "baked in" Q-value
  // already incorporated) and the treatment in Guille's thesis. Differences
  // between the two lead to a few-tens-of-MeV shift in the energy transfer
  // distribution for EM MEC. The shift is done in terms of the binding energy
  // value associated with the original tensor (Eb_ten). Corrections for
  // scaling to a different target are already handled above.
  // - S. Gardiner, 1 July 2020
  bool isEM = interaction.ProcInfo().IsEM();
  if ( isEM ) Delta_Q_value -= 2. * Eb_ten;

  return Delta_Q_value ;
}
//_________________________________________________________________________
void MartiniEricsonChanfrayMarteauMECPXSec2024::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MartiniEricsonChanfrayMarteauMECPXSec2024::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//_________________________________________________________________________
void MartiniEricsonChanfrayMarteauMECPXSec2024::LoadConfig(void)
{
  bool good_config = true ;
  // Cross section scaling factor
  GetParamDef("MEC-XSecScale", fXSecScale, 1.) ;

  fHadronTensorModel = dynamic_cast<const HadronTensorModelI*> ( this->SubAlg("HadronTensorAlg") );
  if( !fHadronTensorModel ) {
    good_config = false ;
    LOG("MartiniEricsonChanfrayMarteauMECPXSec2024", pERROR) << "The required HadronTensorAlg does not exist. AlgoID is : " << SubAlg("HadronTensorAlg")->Id();
  }

  fXSecIntegrator = dynamic_cast<const XSecIntegratorI*> (this->SubAlg("NumericalIntegrationAlg"));
  if( !fXSecIntegrator ) {
    good_config = false ;
    LOG("MartiniEricsonChanfrayMarteauMECPXSec2024", pERROR) << "The required NumericalIntegrationAlg does not exist. AlgId is : " << SubAlg("NumericalIntegrationAlg")->Id() ;
  }

  //Fermi momentum tables for scaling
  this->GetParam( "FermiMomentumTable", fKFTable);

  //binding energy lookups for scaling
  this->GetParam( "RFG-NucRemovalE@Pdg=1000020040", fEbHe );
  this->GetParam( "RFG-NucRemovalE@Pdg=1000030060", fEbLi );
  this->GetParam( "RFG-NucRemovalE@Pdg=1000060120", fEbC  );
  this->GetParam( "RFG-NucRemovalE@Pdg=1000080160", fEbO  );
  this->GetParam( "RFG-NucRemovalE@Pdg=1000120240", fEbMg );
  this->GetParam( "RFG-NucRemovalE@Pdg=1000180400", fEbAr );
  this->GetParam( "RFG-NucRemovalE@Pdg=1000200400", fEbCa );
  this->GetParam( "RFG-NucRemovalE@Pdg=1000260560", fEbFe );
  this->GetParam( "RFG-NucRemovalE@Pdg=1000280580", fEbNi );
  this->GetParam( "RFG-NucRemovalE@Pdg=1000501190", fEbSn );
  this->GetParam( "RFG-NucRemovalE@Pdg=1000791970", fEbAu );
  this->GetParam( "RFG-NucRemovalE@Pdg=1000822080", fEbPb );

  // Read optional MECScaleVsW:
  fMECScaleAlg = nullptr;
  if( GetConfig().Exists("MECScaleAlg") ) {
    fMECScaleAlg = dynamic_cast<const XSecScaleI *> ( this->SubAlg("MECScaleAlg") );
    if( !fMECScaleAlg ) {
      good_config = false ;
      LOG("MartiniEricsonChanfrayMarteauMECPXSec2024", pERROR) << "The required MECScaleAlg cannot be casted. AlgID is : " << SubAlg("MECScaleAlg")->Id() ;
    }
  }

  // Read optional QvalueShifter:
  fQvalueShifter = nullptr;
  if( GetConfig().Exists("QvalueShifterAlg") ) {
    fQvalueShifter = dynamic_cast<const QvalueShifter *> ( this->SubAlg("QvalueShifterAlg") );
    if( !fQvalueShifter ) {
      good_config = false ;
      LOG("MartiniEricsonChanfrayMarteauMECPXSec2024", pERROR) << "The required QvalueShifterAlg does not exist. AlgId is : " << SubAlg("QvalueShifterAlg")->Id() ;
    }
  }

  if( ! good_config ) {
    LOG("MartiniEricsonChanfrayMarteauMECPXSec2024", pERROR) << "Configuration has failed.";
    exit(78) ;
  }

}
