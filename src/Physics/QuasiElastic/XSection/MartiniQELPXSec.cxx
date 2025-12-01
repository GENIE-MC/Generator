//_________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE
 For the class documentation see the corresponding header file.
*/
//_________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/HadronTensors/LabFrameHadronTensorI.h"
#include "Physics/HadronTensors/MartiniQELHadronTensorModel.h"
#include "Physics/QuasiElastic/XSection/MartiniQELPXSec.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
using namespace genie;

//_________________________________________________________________________
MartiniQELPXSec::MartiniQELPXSec() : XSecAlgorithmI("genie::MartiniQELPXSec")
{
}
//_________________________________________________________________________
MartiniQELPXSec::MartiniQELPXSec(string config)
  : XSecAlgorithmI("genie::MartiniQELPXSec", config)
{
}
//_________________________________________________________________________
MartiniQELPXSec::~MartiniQELPXSec()
{
}
//_________________________________________________________________________
double MartiniQELPXSec::XSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{
  if ( !this->ValidProcess(interaction) ) return 0.;

  // Get the hadron tensor for the selected nuclide. Check the probe PDG code
  // to know whether to use the tensor for CC neutrino scattering or for
  // electron scattering
  int target_pdg = interaction->InitState().Tgt().Pdg();
  int probe_pdg = interaction->InitState().ProbePdg();
  int tensor_pdg = target_pdg;
  int A_request = pdg::IonPdgCodeToA(target_pdg);
  bool need_to_scale = false;

  HadronTensorType_t tensor_type = kHT_Undefined;
  if ( pdg::IsNeutrino(probe_pdg) || pdg::IsAntiNeutrino(probe_pdg) ) {
    tensor_type = kHT_QE_Full;
  }
  else {
    // If the probe is not a neutrino, assume that it's an electron
    tensor_type = kHT_QE_EM;
  }

  double Eb_tgt=0;
  double Eb_ten=0;

  if ( A_request <= 4 ) {
    // Use carbon tensor for very light nuclei. This is not ideal . . .
    tensor_pdg = kPdgTgtC12;
    Eb_tgt=fEbHe; Eb_ten=fEbC;
  }
  else if (A_request < 9) {
    tensor_pdg = kPdgTgtC12;
    Eb_tgt=fEbLi; Eb_ten=fEbC;
  }
  else if (A_request >= 9 && A_request < 15) {
    tensor_pdg = kPdgTgtC12;
    Eb_tgt=fEbC; Eb_ten=fEbC;
  }
  else if(A_request >= 15 && A_request < 22) {
    tensor_pdg = kPdgTgtO16;
    Eb_tgt=fEbO; Eb_ten=fEbO;
  }
  else if(A_request >= 22 && A_request < 40) {
    tensor_pdg = kPdgTgtC12;
    Eb_tgt=fEbMg; Eb_ten=fEbO;
  }
  else if(A_request >= 40 && A_request < 56) {
    tensor_pdg = kPdgTgtCa40;
    Eb_tgt=fEbCa; Eb_ten=fEbCa;
  }
  else if(A_request >= 56 && A_request < 119) {
    tensor_pdg = kPdgTgtCa40;
    Eb_tgt=fEbFe; Eb_ten=fEbCa;
  }
  else if(A_request >= 119 && A_request < 206) {
    tensor_pdg = kPdgTgtCa40;
    Eb_tgt=fEbSn; Eb_ten=fEbCa;
  }
  else if(A_request >= 206) {
    tensor_pdg = kPdgTgtCa40;
    Eb_tgt=fEbPb; Eb_ten=fEbCa;
  }

  if (tensor_pdg != target_pdg) need_to_scale = true;
  // The Martini-1p1h hadron tensors are defined using the same conventions
  // as the Valencia MEC (and Martini-MEC) model, so we can use the same sort of tensor
  // object to describe them.
  const LabFrameHadronTensorI* tensor
    = dynamic_cast<const LabFrameHadronTensorI*>( fHadronTensorModel->GetTensor(tensor_pdg,
    tensor_type) );

  // If retrieving the tensor failed, complain and return zero
  if ( !tensor ) {
    LOG("MartiniQE", pWARN) << "Failed to load a hadronic tensor for the"
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

  // SD: The Q-Value essentially corrects q0 to account for nuclear
  // binding energy in the Valencia model but this effect is already
  // in Guille's tensors so I'll set it to 0.
  // However, if I want to scale I need to account for the altered
  // binding energy. To first order I can use the Delta_Q_value for this
  double Delta_Q_value = Eb_tgt-Eb_ten;

  // Apply Qvalue relative shift if needed:
  if( fQvalueShifter ) {
    // We have the option to add an additional shift on top of the binding energy correction
    // The QvalueShifter, is a relative shift to the Q_value. 
    // The Q_value was already taken into account in the hadron tensor. Here we recalculate it
    // to get the right absolute shift. 
    double tensor_Q_value = genie::utils::mec::Qvalue(tensor_pdg,probe_pdg);
    double total_Q_value = tensor_Q_value + Delta_Q_value ; 
    double Q_value_shift = total_Q_value * fQvalueShifter -> Shift( interaction->InitState().Tgt() ) ; 
    Delta_Q_value += Q_value_shift ;
  }

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

  // Compute the cross section using the hadron tensor
  double xsec = tensor->dSigma_dT_dCosTheta_rosenbluth(interaction, Delta_Q_value);
  LOG("MartiniQE", pDEBUG) << "XSec in cm2 / neutron is  " << xsec/(units::cm2);

  // Currently the Martini QE hadron tensors are given per active nucleon, but
  // the calculation above assumes they are per atom. Need to adjust for this.
  const ProcessInfo& proc_info = interaction->ProcInfo();

  // Neutron, proton, and mass numbers of the target
  const Target& tgt = interaction->InitState().Tgt();
  if ( proc_info.IsWeakCC() ) {
    if ( pdg::IsNeutrino(probe_pdg) ) xsec *= tgt.N();
    else if ( pdg::IsAntiNeutrino(probe_pdg) ) xsec *= tgt.Z();
    else {
      // We should never get here if ValidProcess() is working correctly
      LOG("MartiniQE", pERROR) << "Unrecognized probe " << probe_pdg
        << " encountered for a WeakCC process";
      xsec = 0.;
    }
  }
  else if ( proc_info.IsEM() || proc_info.IsWeakNC() ) {
    // For EM processes, scale by the number of nucleons of the same type
    // as the struck one. This ensures the correct ratio of initial-state
    // p vs. n when making splines. The nuclear cross section is obtained
    // by scaling by A/2 for an isoscalar target, so we can get the right
    // behavior for all targets by scaling by Z/2 or N/2 as appropriate.
    // Do the same for NC.
    int hit_nuc_pdg = tgt.HitNucPdg();
    if ( pdg::IsProton(hit_nuc_pdg) ) xsec *= tgt.Z() / 2.;
    else if ( pdg::IsNeutron(hit_nuc_pdg) ) xsec *= tgt.N() / 2.;
    // We should never get here if ValidProcess() is working correctly
    else return 0.;
  }
  else {
    // We should never get here if ValidProcess() is working correctly
    LOG("MartiniQE", pERROR) << "Unrecognized process " << proc_info.AsString()
      << " encountered in MartiniQELPXSec::XSec()";
    xsec = 0.;
  }

  LOG("MartiniQE", pDEBUG) << "XSec in cm2 / atom is  " << xsec / units::cm2;

  // This scaling should be okay-ish for the total xsec, but it misses
  // the energy shift. To get this we should really just build releveant
  // hadron tensors but there may be some ways to approximate it.
  // For more details see Guille's thesis: https://idus.us.es/xmlui/handle/11441/74826
  if ( need_to_scale ) {
    FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
    const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
    double KF_tgt = kft->FindClosestKF(target_pdg, kPdgProton);
    double KF_ten = kft->FindClosestKF(tensor_pdg, kPdgProton);
    LOG("MartiniQE", pDEBUG) << "KF_tgt = " << KF_tgt;
    LOG("MartiniQE", pDEBUG) << "KF_ten = " << KF_ten;
    double scaleFact = (KF_ten/KF_tgt); // A-scaling already applied in section above
    xsec *= scaleFact;
  }

  // Apply given overall scaling factor
  xsec *= fXSecScale;

  if ( kps != kPSTlctl ) {
    LOG("MartiniQE", pWARN)
      << "Doesn't support transformation from "
      << KinePhaseSpace::AsString(kPSTlctl) << " to "
      << KinePhaseSpace::AsString(kps);
    xsec = 0.;
  }

  return xsec;
}
//_________________________________________________________________________
double MartiniQELPXSec::Integral(const Interaction* interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this, interaction);
  return xsec;
}
//_________________________________________________________________________
bool MartiniQELPXSec::ValidProcess(const Interaction* interaction) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if ( !proc_info.IsQuasiElastic() ) return false;

  // The Martini calculation is only appropriate for complex nuclear targets,
  // not free nucleons.
  if ( !init_state.Tgt().IsNucleus() ) return false;

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  bool isP   = pdg::IsProton(nuc);
  bool isN   = pdg::IsNeutron(nuc);
  bool isnu  = pdg::IsNeutrino(nu);
  bool isnub = pdg::IsAntiNeutrino(nu);
  bool is_chgl = pdg::IsChargedLepton(nu);

  bool prcok = ( proc_info.IsWeakCC() && ((isP && isnub) || (isN && isnu)) )
    || ( proc_info.IsEM() && is_chgl && (isP || isN) );
  if ( !prcok ) return false;

  return true;
}
//_________________________________________________________________________
void MartiniQELPXSec::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MartiniQELPXSec::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//_________________________________________________________________________
void MartiniQELPXSec::LoadConfig(void)
{
  bool good_config = true ;

  // Cross section scaling factor
  GetParam( "QEL-CC-XSecScale", fXSecScale ) ;

  fHadronTensorModel = dynamic_cast< const MartiniQELHadronTensorModel* >(
    this->SubAlg("HadronTensorAlg") );
  assert( fHadronTensorModel );

  // Load XSec Integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *>(
    this->SubAlg("XSec-Integrator") );
  assert( fXSecIntegrator );

  // Fermi momentum tables for scaling
  this->GetParam( "FermiMomentumTable", fKFTable);

  // Binding energy lookups for scaling
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

  // Read optional QvalueShifter:
  // Read optional QvalueShifter:                                                                                   
  fQvalueShifter = nullptr;                                                                                        
  if( GetConfig().Exists("QvalueShifterAlg") ) {            
    
    fQvalueShifter = dynamic_cast<const QvalueShifter *> ( this->SubAlg("QvalueShifterAlg") );      

    if( !fQvalueShifter ) {                                                      

      good_config = false ;                                                    
                                  
      LOG("MartiniQE", pERROR) << "The required QvalueShifterAlg is not valid. AlgID is : " 
			      << SubAlg("QvalueShifterAlg")->Id() ;    
    }                                                                                                               
  }  // if there is a requested QvalueShifteralgo                                                                                                                
  if( ! good_config ) {
    LOG("MartiniQE", pERROR) << "Configuration has failed.";
    exit(78) ;
  }
 

}
