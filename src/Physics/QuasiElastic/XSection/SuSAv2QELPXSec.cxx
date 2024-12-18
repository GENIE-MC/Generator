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
#include "Physics/HadronTensors/SuSAv2QELHadronTensorModel.h"
#include "Physics/QuasiElastic/XSection/SuSAv2QELPXSec.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Conventions/Constants.h"

#include <TMath.h>


using namespace genie;
using namespace std::complex_literals;

//_________________________________________________________________________
SuSAv2QELPXSec::SuSAv2QELPXSec() : XSecAlgorithmI("genie::SuSAv2QELPXSec")
{
}
//_________________________________________________________________________
SuSAv2QELPXSec::SuSAv2QELPXSec(string config)
	: XSecAlgorithmI("genie::SuSAv2QELPXSec", config)
{
}
//_________________________________________________________________________
SuSAv2QELPXSec::~SuSAv2QELPXSec()
{
}
//_________________________________________________________________________
double SuSAv2QELPXSec::XSec(const Interaction* interaction,
		KinePhaseSpace_t kps) const
{
	if ( !this->ValidProcess(interaction) ) return 0.;

	// Check that the input kinematical point is within the range
	// in which hadron tensors are known (for chosen target)
	double Ev    = interaction->InitState().ProbeE(kRfLab);
	double Tl    = interaction->Kine().GetKV(kKVTl);
	double costl = interaction->Kine().GetKV(kKVctl);
	double ml    = interaction->FSPrimLepton()->Mass();
	double Q0    = 0.;
	double Q3    = 0.;

	genie::utils::mec::Getq0q3FromTlCostl(Tl, costl, Ev, ml, Q0, Q3);

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


	// ******************************
	// Now choose which tesor to use
	// ******************************

	// Get the hadron tensor for the selected nuclide. Check the probe PDG code
	// to know whether to use the tensor for CC neutrino scattering or for
	// electron scattering
	int target_pdg = interaction->InitState().Tgt().Pdg();
	int probe_pdg = interaction->InitState().ProbePdg();
	// get the hit nucleon pdg code
	int hit_nuc_pdg = interaction->InitState().Tgt().HitNucPdg();

	int tensor_pdg_susa = target_pdg;
	int tensor_pdg_crpa = target_pdg;
	int A_request = pdg::IonPdgCodeToA(target_pdg);
	int Z_request = pdg::IonPdgCodeToZ(target_pdg);
	bool need_to_scale_susa = false;
	bool need_to_scale_crpa = false;


	// This gets a bit messy as the different models have different
	// targets available and currently only SuSA does EM

	HadronTensorType_t tensor_type_susa = kHT_Undefined;
	HadronTensorType_t tensor_type_crpa = kHT_Undefined;
	HadronTensorType_t tensor_type_blen = kHT_Undefined;

	if ( pdg::IsNeutrino(probe_pdg) ) {
		tensor_type_susa = kHT_QE_Full;
		tensor_type_blen = kHT_QE_SuSABlend;
		// CRPA/HF tensors having q0 dependent binning, so are split
		// CRPA
		if (modelConfig == kMd_CRPA || modelConfig == kMd_CRPASuSAv2Hybrid){
			if(Q0<0.060) tensor_type_crpa = kHT_QE_CRPA_Low;
			else if(Q0<0.150) tensor_type_crpa = kHT_QE_CRPA_Medium;
			else tensor_type_crpa = kHT_QE_CRPA_High;
		}
		// Hartree-Fock
		if (modelConfig == kMd_HF || modelConfig == kMd_HFSuSAv2Hybrid){
			if(Q0<0.060) tensor_type_crpa = kHT_QE_HF_Low;
			else if(Q0<0.150) tensor_type_crpa = kHT_QE_HF_Medium;
			else tensor_type_crpa = kHT_QE_HF_High;
		}
		if (modelConfig == kMd_CRPAPW || modelConfig == kMd_CRPAPWSuSAv2Hybrid){
			if(Q0<0.060) tensor_type_crpa = kHT_QE_CRPAPW_Low;
			else if(Q0<0.150) tensor_type_crpa = kHT_QE_CRPAPW_Medium;
			else tensor_type_crpa = kHT_QE_CRPAPW_High;
		}
		// Hartree-Fock
		if (modelConfig == kMd_HFPW || modelConfig == kMd_HFPWSuSAv2Hybrid){
			if(Q0<0.060) tensor_type_crpa = kHT_QE_HFPW_Low;
			else if(Q0<0.150) tensor_type_crpa = kHT_QE_HFPW_Medium;
			else tensor_type_crpa = kHT_QE_HFPW_High;
		}
	}
	else if ( pdg::IsAntiNeutrino(probe_pdg) ){
		// SuSA implementation doesn't accoutn for asymmetry between protons
		// and neutrons. In general this is a small effect.
		tensor_type_susa = kHT_QE_Full;
		//For the blending case, Ar40 is treated specially:
		if(A_request == 40 && Z_request == 18){
			tensor_type_blen = kHT_QE_SuSABlend_anu;
		}
		else tensor_type_blen = kHT_QE_SuSABlend;
		// CRPA tensors having q0 dependent binning, so are split:
		//CRPA
		if (modelConfig == kMd_CRPA || modelConfig == kMd_CRPASuSAv2Hybrid){
			if(Q0<0.060) tensor_type_crpa = kHT_QE_CRPA_anu_Low;
			else if(Q0<0.150) tensor_type_crpa = kHT_QE_CRPA_anu_Medium;
			else tensor_type_crpa = kHT_QE_CRPA_anu_High;
		}
		// Hartree-Fock
		if (modelConfig == kMd_HF || modelConfig == kMd_HFSuSAv2Hybrid){
			if(Q0<0.060) tensor_type_crpa = kHT_QE_HF_anu_Low;
			else if(Q0<0.150) tensor_type_crpa = kHT_QE_HF_anu_Medium;
			else tensor_type_crpa = kHT_QE_HF_anu_High;
		}
		if (modelConfig == kMd_CRPAPW || modelConfig == kMd_CRPAPWSuSAv2Hybrid){
			if(Q0<0.060) tensor_type_crpa = kHT_QE_CRPAPW_anu_Low;
			else if(Q0<0.150) tensor_type_crpa = kHT_QE_CRPAPW_anu_Medium;
			else tensor_type_crpa = kHT_QE_CRPAPW_anu_High;
		}
		// Hartree-Fock
		if (modelConfig == kMd_HFPW || modelConfig == kMd_HFPWSuSAv2Hybrid){
			if(Q0<0.060) tensor_type_crpa = kHT_QE_HFPW_anu_Low;
			else if(Q0<0.150) tensor_type_crpa = kHT_QE_HFPW_anu_Medium;
			else tensor_type_crpa = kHT_QE_HFPW_anu_High;
		}
	}
	else {
		// If the probe is not a neutrino, assume that it's an electron
		// Currently only avaialble for SuSA. CRPA coming soon(ish)!
		if ( pdg::IsProton(hit_nuc_pdg) ) {
			tensor_type_susa = kHT_QE_EM_proton;
		}
		else if( pdg::IsNeutron(hit_nuc_pdg) ) {
			tensor_type_susa = kHT_QE_EM_neutron;
		}
		else {
			tensor_type_susa = kHT_QE_EM;  	// default
		}

	}

	double Eb_tgt=0;
	double Eb_ten_susa=0;
	double Eb_ten_crpa=0;

	if ( A_request <= 4 ) {
		// Use carbon tensor for very light nuclei. This is not ideal . . .
		Eb_tgt = fEbHe;
		tensor_pdg_susa = kPdgTgtC12;
		tensor_pdg_crpa = kPdgTgtC12;
		Eb_ten_susa = fEbC;
		Eb_ten_crpa = fEbC;
	}
	else if (A_request < 9) {
		Eb_tgt=fEbLi;
		tensor_pdg_susa = kPdgTgtC12;
		tensor_pdg_crpa = kPdgTgtC12;
		Eb_ten_susa = fEbC;
		Eb_ten_crpa = fEbC;
	}
	else if (A_request >= 9 && A_request < 15) {
		Eb_tgt=fEbC;
		tensor_pdg_susa = kPdgTgtC12;
		tensor_pdg_crpa = kPdgTgtC12;
		Eb_ten_susa = fEbC;
		Eb_ten_crpa = fEbC;
	}
	else if(A_request >= 15 && A_request < 22) {
		Eb_tgt=fEbO;
		tensor_pdg_susa = kPdgTgtC12;
		tensor_pdg_crpa = kPdgTgtO16;
		Eb_ten_susa = fEbC;
		Eb_ten_crpa = fEbO;
	}
	else if(A_request == 40 && Z_request == 18) {
		// Treat the common non-isoscalar case specially
		Eb_tgt=fEbAr;
		tensor_pdg_susa = kPdgTgtC12;
		tensor_pdg_crpa = 1000180400;
		Eb_ten_susa = fEbC;
		Eb_ten_crpa = fEbAr;
	}
	else if(A_request >= 22 && A_request < 40) {
		Eb_tgt=fEbMg;
		tensor_pdg_susa = kPdgTgtC12;
		tensor_pdg_crpa = kPdgTgtO16;
		Eb_ten_susa = fEbC;
		Eb_ten_crpa = fEbO;
	}
	else if(A_request >= 40 && A_request < 56) {
		Eb_tgt=fEbAr;
		tensor_pdg_susa = kPdgTgtC12;
		tensor_pdg_crpa = kPdgTgtO16;
		Eb_ten_susa = fEbC;
		Eb_ten_crpa = fEbO;
	}
	else if(A_request >= 56 && A_request < 119) {
		Eb_tgt=fEbFe;
		tensor_pdg_susa = kPdgTgtC12;
		tensor_pdg_crpa = kPdgTgtO16;
		Eb_ten_susa = fEbC;
		Eb_ten_crpa = fEbO;
	}
	else if(A_request >= 119 && A_request < 206) {
		Eb_tgt=fEbSn;
		tensor_pdg_susa = kPdgTgtC12;
		tensor_pdg_crpa = kPdgTgtO16;
		Eb_ten_susa = fEbC;
		Eb_ten_crpa = fEbO;
	}
	else if(A_request >= 206) {
		Eb_tgt=fEbPb;
		tensor_pdg_susa = kPdgTgtC12;
		tensor_pdg_crpa = kPdgTgtO16;
		Eb_ten_susa = fEbC;
		Eb_ten_crpa = fEbO;
	}

	if (tensor_pdg_susa != target_pdg) need_to_scale_susa = true;
	if (tensor_pdg_crpa != target_pdg) need_to_scale_crpa = true;

	// Finally we can now get the tensors we need

	const LabFrameHadronTensorI* tensor_susa;
	const LabFrameHadronTensorI* tensor_crpa;
	const LabFrameHadronTensorI* tensor_blen;

	if( modelConfig == kMd_SuSAv2 ){
		tensor_susa = dynamic_cast<const LabFrameHadronTensorI*>
			( fHadronTensorModel->GetTensor (tensor_pdg_susa, tensor_type_susa) );

		if ( !tensor_susa ) {
			LOG("SuSAv2QE", pWARN) << "Failed to load a SuSAv2 hadronic tensor for the"
				" nuclide " << tensor_pdg_susa;
			return 0.;
		}
	}

	if( modelConfig == kMd_CRPASuSAv2Hybrid   || modelConfig == kMd_HFSuSAv2Hybrid ||
			modelConfig == kMd_CRPAPWSuSAv2Hybrid || modelConfig == kMd_HFPWSuSAv2Hybrid ||
			modelConfig == kMd_SuSAv2Blend){
		tensor_blen = dynamic_cast<const LabFrameHadronTensorI*>
			( fHadronTensorModel->GetTensor (tensor_pdg_crpa, tensor_type_blen) );

		// If retrieving the tensor failed, complain and return zero
		if ( !tensor_blen ) {
			LOG("SuSAv2QE", pWARN) << "Failed to load a blending SuSAv2 hadronic tensor for the"
				" nuclide " << tensor_pdg_crpa;
			return 0.;
		}
	}

	if( modelConfig != kMd_SuSAv2 && modelConfig != kMd_SuSAv2Blend){

		tensor_crpa = dynamic_cast<const LabFrameHadronTensorI*>
			( fHadronTensorModel->GetTensor (tensor_pdg_crpa, tensor_type_crpa) );

		// If retrieving the tensor failed, complain and return zero
		if ( !tensor_crpa ) {
			LOG("SuSAv2QE", pWARN) << "Failed to load a CRPA or HF hadronic tensor for the"
				" nuclide " << tensor_pdg_crpa;
			return 0.;
		}
	}

	// *****************************
	// Q_value offset calculation
	// *****************************

	// SD: The Q-Value essentially corrects q0 to account for nuclear
	// binding energy in the Valencia model but this effect is already
	// in the SuSAv2 and CRPA/HF tensors so I'll set it to 0.
	// However, if I want to scale I need to account for the altered
	// binding energy. To first order I can use the Q_value for this
	double Delta_Q_value_susa = Eb_tgt-Eb_ten_susa;
	double Delta_Q_value_crpa = Eb_tgt-Eb_ten_crpa;
	double Delta_Q_value_blen = Eb_tgt-Eb_ten_crpa;

	// Apply Qvalue relative shift if needed:
	if( fQvalueShifter ) {
		// We have the option to add an additional shift on top of the binding energy correction
		// The QvalueShifter, is a relative shift to the Q_value.
		// The Q_value was already taken into account in the hadron tensor. Here we recalculate it
		// to get the right absolute shift.
		double tensor_Q_value_susa = genie::utils::mec::Qvalue(tensor_pdg_susa,probe_pdg);
		double total_Q_value_susa = tensor_Q_value_susa + Delta_Q_value_susa ;
		double Q_value_shift_susa = total_Q_value_susa * fQvalueShifter -> Shift( interaction->InitState().Tgt() ) ;

		double tensor_Q_value_crpa = genie::utils::mec::Qvalue(tensor_pdg_crpa,probe_pdg);
		double total_Q_value_crpa = tensor_Q_value_crpa + Delta_Q_value_crpa ;
		double Q_value_shift_crpa = total_Q_value_crpa * fQvalueShifter -> Shift( interaction->InitState().Tgt() ) ;

		Delta_Q_value_susa += Q_value_shift_susa;
		Delta_Q_value_crpa += Q_value_shift_crpa;
		Delta_Q_value_blen += Q_value_shift_crpa;
	}

	// Set the xsec to zero for interactions with q0,q3 outside the requested range

	if( modelConfig == kMd_SuSAv2){
		double Q0min = tensor_susa->q0Min();
		double Q0max = tensor_susa->q0Max();
		double Q3min = tensor_susa->qMagMin();
		double Q3max = tensor_susa->qMagMax();
		if (Q0-Delta_Q_value_susa < Q0min || Q0-Delta_Q_value_susa > Q0max || Q3 < Q3min || Q3 > Q3max) {
			return 0.0;
		}
	}

	else if ( modelConfig == kMd_CRPA   || modelConfig == kMd_HF ||
			modelConfig == kMd_CRPAPW || modelConfig == kMd_HFPW ){
		double Q0min = tensor_crpa->q0Min();
		double Q0max = tensor_crpa->q0Max();
		double Q3min = tensor_crpa->qMagMin();
		double Q3max = tensor_crpa->qMagMax();
		if (Q0-Delta_Q_value_crpa < Q0min || Q0-Delta_Q_value_crpa > Q0max || Q3 < Q3min || Q3 > Q3max) {
			return 0.0;
		}
	}

	else if ( modelConfig == kMd_SuSAv2Blend){
		double Q0min = tensor_blen->q0Min();
		double Q0max = tensor_blen->q0Max();
		double Q3min = tensor_blen->qMagMin();
		double Q3max = tensor_blen->qMagMax();
		if (Q0-Delta_Q_value_blen < Q0min || Q0-Delta_Q_value_blen > Q0max || Q3 < Q3min || Q3 > Q3max) {
			return 0.0;
		}
	}

	else{ // hybrid (blending) cases. Low kinematics handled by CRPA/HF, high kinematics by blended SuSA
		double Q0min = tensor_crpa->q0Min();
		double Q0max = tensor_blen->q0Max();
		double Q3min = tensor_crpa->qMagMin();
		double Q3max = tensor_blen->qMagMax();
		if (Q0-Delta_Q_value_crpa < Q0min || Q0-Delta_Q_value_blen > Q0max || Q3 < Q3min || Q3 > Q3max) {
			return 0.0;
		}
	}

	// ******************************
	// Actual xsec calculation
	// ******************************

	double xsec_susa = 0;
	double xsec_crpa = 0;
	double xsec_blen = 0;

	if( modelConfig == kMd_SuSAv2 ){
		// Compute the cross section using the hadron tensor
		xsec_susa = tensor_susa->dSigma_dT_dCosTheta_rosenbluth(interaction, Delta_Q_value_susa);
		LOG("SuSAv2QE", pDEBUG) << "SuSAv2 XSec in cm2 / neutron is  " << xsec_susa/(units::cm2);
		xsec_susa = XSecScaling(xsec_susa, interaction, target_pdg, tensor_pdg_susa, need_to_scale_susa);
	}


	if( modelConfig == kMd_CRPASuSAv2Hybrid   || modelConfig == kMd_HFSuSAv2Hybrid ||
			modelConfig == kMd_CRPAPWSuSAv2Hybrid || modelConfig == kMd_HFPWSuSAv2Hybrid ||
			modelConfig == kMd_SuSAv2Blend){
		// Compute the cross section using the hadron tensor
		xsec_blen = tensor_blen->dSigma_dT_dCosTheta_rosenbluth(interaction, Delta_Q_value_blen);
		LOG("SuSAv2QE", pDEBUG) << "SuSAv2 (blending) XSec in cm2 / atom is  " << xsec_blen/(units::cm2);
		// The blended SuSAv2 calculation already gives the xsec per atom
		// For the A-scaling below to make sense we need to transform them to per active nucleon
		int A_tensor = pdg::IonPdgCodeToA(tensor_pdg_crpa);
		int Z_tensor = pdg::IonPdgCodeToZ(tensor_pdg_crpa);
		int N_tensor = A_tensor-Z_tensor;

		if ( pdg::IsNeutrino(probe_pdg) ) xsec_blen *= 1.0/N_tensor;
		else if ( pdg::IsAntiNeutrino(probe_pdg) ) xsec_blen *= 1.0/Z_tensor;

		xsec_blen = XSecScaling(xsec_blen, interaction, target_pdg, tensor_pdg_crpa, need_to_scale_crpa);
	}

	if( modelConfig != kMd_SuSAv2 && modelConfig != kMd_SuSAv2Blend){
		// Compute the cross section using the hadron tensor
		xsec_crpa = tensor_crpa->dSigma_dT_dCosTheta_rosenbluth(interaction, Delta_Q_value_crpa);
		LOG("SuSAv2QE", pDEBUG) << "CRPA or HF XSec in cm2 / atom is  " << xsec_crpa/(units::cm2);
		// The CRPA calculation already gives the xsec per atom
		// For the A-scaling below to make sense we need to transform them to per active nucleon
		int A_tensor = pdg::IonPdgCodeToA(tensor_pdg_crpa);
		int Z_tensor = pdg::IonPdgCodeToZ(tensor_pdg_crpa);
		int N_tensor = A_tensor-Z_tensor;

		if ( pdg::IsNeutrino(probe_pdg) ) xsec_crpa *= 1.0/N_tensor;
		else if ( pdg::IsAntiNeutrino(probe_pdg) ) xsec_crpa *= 1.0/Z_tensor;

		// TODO: When we add EM for CRPA need to add how this case it dealt with here

		xsec_crpa = XSecScaling(xsec_crpa, interaction, target_pdg, tensor_pdg_crpa, need_to_scale_crpa);

	}

	// Apply blending if needed
	double xsec = 0;

	if( modelConfig == kMd_SuSAv2 )      xsec = xsec_susa;
	if( modelConfig == kMd_SuSAv2Blend ) xsec = xsec_blen;
	if( modelConfig == kMd_CRPA   || modelConfig == kMd_HF ||
			modelConfig == kMd_CRPAPW || modelConfig == kMd_HFPW ) xsec = xsec_crpa;
	else if( modelConfig == kMd_CRPASuSAv2Hybrid   ||
			modelConfig == kMd_HFSuSAv2Hybrid     ||
			modelConfig == kMd_CRPAPWSuSAv2Hybrid ||
			modelConfig == kMd_HFPWSuSAv2Hybrid ){  // blending cases
		if(blendMode == 1) // Linear blending in q0
			if      (Q0 < q0BlendStart)  xsec = xsec_crpa;
			else if (Q0 > q0BlendEnd)    xsec = xsec_blen;
			else{
				double SuSAFrac = (Q0 - q0BlendStart) / (q0BlendEnd - q0BlendStart);
				double CRPAFrac = 1 - SuSAFrac;
				xsec = SuSAFrac*xsec_blen + CRPAFrac*xsec_crpa;
				LOG("SuSAv2QE", pDEBUG) << "Q0 is  " << Q0;
				LOG("SuSAv2QE", pDEBUG) << "SuSAFrac is  " << SuSAFrac;
				LOG("SuSAv2QE", pDEBUG) << "CRPAFrac is  " << CRPAFrac;
				LOG("SuSAv2QE", pDEBUG) << "xsec is  " << xsec;
			}
		else if(blendMode == 2){ // Exp blending in q (from Alexis)
			double phi_q = (genie::constants::kPi / 2.) * (1 - 1./(1+std::exp( (Q3 - qBlendRef)/qBlendDel)) );
			xsec = TMath::Sin(phi_q)*TMath::Sin(phi_q)*xsec_blen + TMath::Cos(phi_q)*TMath::Cos(phi_q)*xsec_crpa;
			LOG("SuSAv2QE", pDEBUG) << "Q3 is  " << Q3;
			LOG("SuSAv2QE", pDEBUG) << "SuSAFrac is  " << TMath::Sin(phi_q)*TMath::Sin(phi_q);
			LOG("SuSAv2QE", pDEBUG) << "CRPAFrac is  " << TMath::Cos(phi_q)*TMath::Cos(phi_q);
			LOG("SuSAv2QE", pDEBUG) << "xsec is  " << xsec;
		}
	}

	// Apply given overall scaling factor
	double xsec_scale = 1 ;
	if( interaction->ProcInfo().IsWeakCC() ) xsec_scale = fXSecCCScale;
	else if( interaction->ProcInfo().IsWeakNC() ) xsec_scale = fXSecNCScale;
	else if( interaction->ProcInfo().IsEM() ) xsec_scale = fXSecEMScale;

	xsec *= xsec_scale ;

    if ( kps != kPSTlctl ) 
    {
      // Compute the appropriate Jacobian for transformation to the requested phase space
      double J = utils::kinematics::Jacobian(interaction, kPSTlctl, kps);
      xsec *= J;
    }

	return xsec;
}

//_________________________________________________________________________
double SuSAv2QELPXSec::XSecScaling(double xsec, const Interaction* interaction, int target_pdg, int tensor_pdg, bool need_to_scale) const
{
	// The xsecs need to be given per active nucleon, but the calculations above are per atom.
	// We also need to A-scale anyhow if the target nucleus is not exactly the one we have the tensor for.
	// We adjust for this bellow.

	const ProcessInfo& proc_info = interaction->ProcInfo();

	// Neutron, proton, and mass numbers of the target
	const Target& tgt = interaction->InitState().Tgt();

	int probe_pdg = interaction->InitState().ProbePdg();

	if ( proc_info.IsWeakCC() ) {
		if ( pdg::IsNeutrino(probe_pdg) ) xsec *= tgt.N();
		else if ( pdg::IsAntiNeutrino(probe_pdg) ) xsec *= tgt.Z();
		else {
			// We should never get here if ValidProcess() is working correctly
			LOG("SuSAv2QE", pERROR) << "Unrecognized probe " << probe_pdg
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
		// Do the same for NC. TODO: double-check that this is the right
		// thing to do when we SuSAv2 NC hadronic tensors are added to GENIE.
		int hit_nuc_pdg = tgt.HitNucPdg();
		if ( pdg::IsProton(hit_nuc_pdg) ) xsec *= tgt.Z() / 2.;
		else if ( pdg::IsNeutron(hit_nuc_pdg) ) xsec *= tgt.N() / 2.;
		// We should never get here if ValidProcess() is working correctly
		else return 0.;
	}
	else {
		// We should never get here if ValidProcess() is working correctly
		LOG("SuSAv2QE", pERROR) << "Unrecognized process " << proc_info.AsString()
			<< " encountered in SuSAv2QELPXSec::XSec()";
		xsec = 0.;
	}

	LOG("SuSAv2QE", pDEBUG) << "XSec in cm2 / atom is  " << xsec / units::cm2;

	// This scaling should be okay-ish for the total xsec, but it misses
	// the energy shift. To get this we should really just build releveant
	// hadron tensors but there may be some ways to approximate it.
	// For more details see Guille's thesis: https://idus.us.es/xmlui/handle/11441/74826

	// We already did some of this when we apply the Q-value shift. We can do a little
	// better by tuning the A-scaling as below, following the SuperScaling ansatz
	if ( need_to_scale ) {
		FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
		const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
		double KF_tgt = kft->FindClosestKF(target_pdg, kPdgProton);
		double KF_ten = kft->FindClosestKF(tensor_pdg, kPdgProton);
		LOG("SuSAv2QE", pDEBUG) << "KF_tgt = " << KF_tgt;
		LOG("SuSAv2QE", pDEBUG) << "KF_ten = " << KF_ten;
		double scaleFact = (KF_ten/KF_tgt); // A-scaling already applied in section above
		xsec *= scaleFact;
	}

	return xsec;
}

//_________________________________________________________________________
double SuSAv2QELPXSec::Integral(const Interaction* interaction) const
{
	double xsec = fXSecIntegrator->Integrate(this, interaction);
	return xsec;
}
//_________________________________________________________________________
bool SuSAv2QELPXSec::ValidProcess(const Interaction* interaction) const
{
	if ( interaction->TestBit(kISkipProcessChk) ) return true;

	const InitialState & init_state = interaction->InitState();
	const ProcessInfo &  proc_info  = interaction->ProcInfo();

	if ( !proc_info.IsQuasiElastic() ) return false;

	// The calculation is only appropriate for complex nuclear targets,
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
void SuSAv2QELPXSec::Configure(const Registry& config)
{
	Algorithm::Configure(config);
	this->LoadConfig();
}
//____________________________________________________________________________
void SuSAv2QELPXSec::Configure(std::string config)
{
	Algorithm::Configure(config);
	this->LoadConfig();
}
//_________________________________________________________________________
void SuSAv2QELPXSec::LoadConfig(void)
{
	bool good_config = true ;

	// Cross section scaling factor
	GetParam( "QEL-CC-XSecScale", fXSecCCScale ) ;
	GetParam( "QEL-NC-XSecScale", fXSecNCScale ) ;
	GetParam( "QEL-EM-XSecScale", fXSecEMScale ) ;
    
  // Do precise calculation of lepton polarization
  GetParamDef( "PreciseLeptonPol", fIsPreciseLeptonPolarization, false ) ;

	// Cross section model choice
	int modelChoice;
	GetParam( "Model-Config", modelChoice ) ;
	modelConfig = (modelType)modelChoice;

	// Blending parameters
	GetParam( "Blend-Mode", blendMode ) ; // 1 = linear, 2 = exp
	GetParam( "q0-Blend-Start", q0BlendStart ) ; // Used for linear
	GetParam( "q0-Blend-End", q0BlendEnd ) ; // Used for linear
	GetParam( "q-Blend-del", qBlendDel ) ; // Used for exp
	GetParam( "q-Blend-ref", qBlendRef ) ; // Used for exp

	fHadronTensorModel = dynamic_cast< const SuSAv2QELHadronTensorModel* >(
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
	fQvalueShifter = nullptr;
	if( GetConfig().Exists("QvalueShifterAlg") ) {

		fQvalueShifter = dynamic_cast<const QvalueShifter *> ( this->SubAlg("QvalueShifterAlg") );

		if( !fQvalueShifter ) {

			good_config = false ;

			LOG("SuSAv2QE", pERROR) << "The required QvalueShifterAlg is not valid. AlgID is : "
				<< SubAlg("QvalueShifterAlg")->Id() ;
		}
	}  // if there is a requested QvalueShifteralgo
	if( ! good_config ) {
		LOG("SuSAv2QE", pERROR) << "Configuration has failed.";
		exit(78) ;
	}

}
const TVector3 & SuSAv2QELPXSec::FinalLeptonPolarization (const Interaction* interaction) const
{    
    const ProcessInfo& proc_info = interaction->ProcInfo();
    if ( proc_info.IsEM() ) 
    {
        LOG("SuSAv2QE", pWARN) << "For EM processes doesn't work yet. Set it to zero.";
        fFinalLeptonPolarization = TVector3(0, 0, 0);
        return fFinalLeptonPolarization;
    }
    if (!fIsPreciseLeptonPolarization || proc_info.IsWeakNC()) return XSecAlgorithmI::FinalLeptonPolarization(interaction);
    const Kinematics&   kinematics = interaction -> Kine();
    const InitialState& init_state = interaction -> InitState();
    const Target& tgt = init_state.Tgt();
  
    bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  
    // HitNucMass() looks up the PDGLibrary (on-shell) value for the initial
    // struck nucleon
    double Mi_onshell = tgt.HitNucMass();
  
    // On-shell mass of final nucleon (from PDGLibrary)
    double Mf = interaction->RecoilNucleon()->Mass();

    // Isoscalar mass of nucleon
    double M = (Mi_onshell + Mf)/2;
    double M2 = M*M;
    
    // Check that the input kinematical point is within the range
    // in which hadron tensors are known (for chosen target)
    double Ev    = init_state.ProbeE(kRfLab);
    double Tl    = kinematics.GetKV(kKVTl);
    double costl = kinematics.GetKV(kKVctl);
    double ml    = interaction->FSPrimLepton()->Mass();
    double Q0    = 0.;
    double Q3    = 0.;

    genie::utils::mec::Getq0q3FromTlCostl(Tl, costl, Ev, ml, Q0, Q3);

    // *** Enforce the global Q^2 cut ***
    double Q2min = genie::controls::kMinQ2Limit; // CC limit

    // Neglect shift due to binding energy. The cut is on the actual
    // value of Q^2, not the effective one to use in the tensor contraction.
    double Q2 = Q3*Q3 - Q0*Q0;
    if ( Q2 < Q2min )
    {
        LOG("SuSAv2QE", pWARN) << "Can't calculate final lepton polarization. Set it to zero.";
        fFinalLeptonPolarization = TVector3(0, 0, 0);
        return fFinalLeptonPolarization;
    }
    
    // Note that GetProbeP4 defaults to returning the probe 4-momentum in the
    // struck nucleon rest frame, so we have to explicitly ask for the lab frame
    // here
    TLorentzVector * tempNeutrino = init_state.GetProbeP4(kRfLab);
    TLorentzVector neutrinoMom = *tempNeutrino;
    delete tempNeutrino;
    TLorentzVector inNucleonMom(*init_state.TgtPtr()->HitNucP4Ptr());
    TLorentzVector inNucleonMomOnShell(inNucleonMom);

    const TLorentzVector leptonMom = kinematics.FSLeptonP4();
    const TLorentzVector outNucleonMom = kinematics.HadSystP4();
    TLorentzVector outNucleonMomOnShell(outNucleonMom);
    
    double inNucleonOnShellEnergy  = TMath::Hypot(M,  inNucleonMomOnShell.P() );
    inNucleonMomOnShell.SetE(inNucleonOnShellEnergy);
    double outNucleonOnShellEnergy = TMath::Hypot(M, outNucleonMomOnShell.P() );
    outNucleonMomOnShell.SetE(outNucleonOnShellEnergy);
    
    // Effective 4-momentum transfer (according to the deForest prescription) for
    // use in computing the hadronic tensor
    TLorentzVector qTildeP4 = outNucleonMomOnShell - inNucleonMomOnShell;
    
    // ******************************
    // Now choose which tesor to use
    // ******************************

    // Get the hadron tensor for the selected nuclide. Check the probe PDG code
    // to use the tensor for CC neutrino scattering
    int target_pdg = interaction->InitState().Tgt().Pdg();
    int probe_pdg = interaction->InitState().ProbePdg();

    int tensor_pdg_susa = target_pdg;
    int tensor_pdg_crpa = target_pdg;
    int A_request = pdg::IonPdgCodeToA(target_pdg);
    int Z_request = pdg::IonPdgCodeToZ(target_pdg);

    // This gets a bit messy as the different models have different
    // targets available and currently only SuSA does EM
    HadronTensorType_t tensor_type_susa = kHT_Undefined;
    HadronTensorType_t tensor_type_crpa = kHT_Undefined;
    HadronTensorType_t tensor_type_blen = kHT_Undefined;

    if ( pdg::IsNeutrino(probe_pdg) ) 
    {
        tensor_type_susa = kHT_QE_Full;
        tensor_type_blen = kHT_QE_SuSABlend;
        // CRPA/HF tensors having q0 dependent binning, so are split
        // CRPA
        if (modelConfig == kMd_CRPA || modelConfig == kMd_CRPASuSAv2Hybrid)
        {
            if(Q0<0.060) tensor_type_crpa = kHT_QE_CRPA_Low;
            else if(Q0<0.150) tensor_type_crpa = kHT_QE_CRPA_Medium;
            else tensor_type_crpa = kHT_QE_CRPA_High;
        }
        // Hartree-Fock
        if (modelConfig == kMd_HF || modelConfig == kMd_HFSuSAv2Hybrid)
        {
            if(Q0<0.060) tensor_type_crpa = kHT_QE_HF_Low;
            else if(Q0<0.150) tensor_type_crpa = kHT_QE_HF_Medium;
            else tensor_type_crpa = kHT_QE_HF_High;
        }
        if (modelConfig == kMd_CRPAPW || modelConfig == kMd_CRPAPWSuSAv2Hybrid)
        {
            if(Q0<0.060) tensor_type_crpa = kHT_QE_CRPAPW_Low;
            else if(Q0<0.150) tensor_type_crpa = kHT_QE_CRPAPW_Medium;
            else tensor_type_crpa = kHT_QE_CRPAPW_High;
        }
        // Hartree-Fock
        if (modelConfig == kMd_HFPW || modelConfig == kMd_HFPWSuSAv2Hybrid)
        {
            if(Q0<0.060) tensor_type_crpa = kHT_QE_HFPW_Low;
            else if(Q0<0.150) tensor_type_crpa = kHT_QE_HFPW_Medium;
            else tensor_type_crpa = kHT_QE_HFPW_High;
        }
    }
    else if ( pdg::IsAntiNeutrino(probe_pdg) )
    {
        // SuSA implementation doesn't accoutn for asymmetry between protons
        // and neutrons. In general this is a small effect.
        tensor_type_susa = kHT_QE_Full;
        //For the blending case, Ar40 is treated specially:
        if(A_request == 40 && Z_request == 18)
        {
            tensor_type_blen = kHT_QE_SuSABlend_anu;
        }
        else tensor_type_blen = kHT_QE_SuSABlend;
        // CRPA tensors having q0 dependent binning, so are split:
        //CRPA
        if (modelConfig == kMd_CRPA || modelConfig == kMd_CRPASuSAv2Hybrid)
        {
            if(Q0<0.060) tensor_type_crpa = kHT_QE_CRPA_anu_Low;
            else if(Q0<0.150) tensor_type_crpa = kHT_QE_CRPA_anu_Medium;
            else tensor_type_crpa = kHT_QE_CRPA_anu_High;
        }
        // Hartree-Fock
        if (modelConfig == kMd_HF || modelConfig == kMd_HFSuSAv2Hybrid)
        {
            if(Q0<0.060) tensor_type_crpa = kHT_QE_HF_anu_Low;
            else if(Q0<0.150) tensor_type_crpa = kHT_QE_HF_anu_Medium;
            else tensor_type_crpa = kHT_QE_HF_anu_High;
        }
        if (modelConfig == kMd_CRPAPW || modelConfig == kMd_CRPAPWSuSAv2Hybrid)
        {
            if(Q0<0.060) tensor_type_crpa = kHT_QE_CRPAPW_anu_Low;
            else if(Q0<0.150) tensor_type_crpa = kHT_QE_CRPAPW_anu_Medium;
            else tensor_type_crpa = kHT_QE_CRPAPW_anu_High;
        }
        // Hartree-Fock
        if (modelConfig == kMd_HFPW || modelConfig == kMd_HFPWSuSAv2Hybrid)
        {
            if(Q0<0.060) tensor_type_crpa = kHT_QE_HFPW_anu_Low;
            else if(Q0<0.150) tensor_type_crpa = kHT_QE_HFPW_anu_Medium;
            else tensor_type_crpa = kHT_QE_HFPW_anu_High;
        }
    }
    else 
    {
        LOG("SuSAv2QE", pWARN) << "For EM processes doesn't work yet. Set it to zero.";
        fFinalLeptonPolarization = TVector3(0, 0, 0);
        return fFinalLeptonPolarization;
    }

    double Eb_tgt=0;
    double Eb_ten_susa=0;
    double Eb_ten_crpa=0;

    if ( A_request <= 4 ) 
    {
        // Use carbon tensor for very light nuclei. This is not ideal . . .
        Eb_tgt = fEbHe;
        tensor_pdg_susa = kPdgTgtC12;
        tensor_pdg_crpa = kPdgTgtC12;
        Eb_ten_susa = fEbC;
        Eb_ten_crpa = fEbC;
    }
    else if (A_request < 9) 
    {
        Eb_tgt=fEbLi;
        tensor_pdg_susa = kPdgTgtC12;
        tensor_pdg_crpa = kPdgTgtC12;
        Eb_ten_susa = fEbC;
        Eb_ten_crpa = fEbC;
    }
    else if (A_request >= 9 && A_request < 15) 
    {
        Eb_tgt=fEbC;
        tensor_pdg_susa = kPdgTgtC12;
        tensor_pdg_crpa = kPdgTgtC12;
        Eb_ten_susa = fEbC;
        Eb_ten_crpa = fEbC;
    }
    else if(A_request >= 15 && A_request < 22) 
    {
        Eb_tgt=fEbO;
        tensor_pdg_susa = kPdgTgtC12;
        tensor_pdg_crpa = kPdgTgtO16;
        Eb_ten_susa = fEbC;
        Eb_ten_crpa = fEbO;
    }
    else if(A_request == 40 && Z_request == 18) 
    {
        // Treat the common non-isoscalar case specially
        Eb_tgt=fEbAr;
        tensor_pdg_susa = kPdgTgtC12;
        tensor_pdg_crpa = 1000180400;
        Eb_ten_susa = fEbC;
        Eb_ten_crpa = fEbAr;
    }
    else if(A_request >= 22 && A_request < 40) 
    {
        Eb_tgt=fEbMg;
        tensor_pdg_susa = kPdgTgtC12;
        tensor_pdg_crpa = kPdgTgtO16;
        Eb_ten_susa = fEbC;
        Eb_ten_crpa = fEbO;
    }
    else if(A_request >= 40 && A_request < 56) {
        Eb_tgt=fEbAr;
        tensor_pdg_susa = kPdgTgtC12;
        tensor_pdg_crpa = kPdgTgtO16;
        Eb_ten_susa = fEbC;
        Eb_ten_crpa = fEbO;
    }
    else if(A_request >= 56 && A_request < 119) 
    {
        Eb_tgt=fEbFe;
        tensor_pdg_susa = kPdgTgtC12;
        tensor_pdg_crpa = kPdgTgtO16;
        Eb_ten_susa = fEbC;
        Eb_ten_crpa = fEbO;
    }
    else if(A_request >= 119 && A_request < 206) 
    {
        Eb_tgt=fEbSn;
        tensor_pdg_susa = kPdgTgtC12;
        tensor_pdg_crpa = kPdgTgtO16;
        Eb_ten_susa = fEbC;
        Eb_ten_crpa = fEbO;
    }
    else if(A_request >= 206) 
    {
        Eb_tgt=fEbPb;
        tensor_pdg_susa = kPdgTgtC12;
        tensor_pdg_crpa = kPdgTgtO16;
        Eb_ten_susa = fEbC;
        Eb_ten_crpa = fEbO;
    }

	// Finally we can now get the tensors we need
	const LabFrameHadronTensorI* tensor_susa;
	const LabFrameHadronTensorI* tensor_crpa;
	const LabFrameHadronTensorI* tensor_blen;

    if( modelConfig == kMd_SuSAv2 )
    {
        tensor_susa = dynamic_cast<const LabFrameHadronTensorI*>
            ( fHadronTensorModel->GetTensor (tensor_pdg_susa, tensor_type_susa) );

        if ( !tensor_susa ) 
        {
            LOG("SuSAv2QE", pWARN) << "Can't calculate final lepton polarization. Set it to zero.";
            fFinalLeptonPolarization = TVector3(0, 0, 0);
            return fFinalLeptonPolarization;
        }
    }

    if( modelConfig == kMd_CRPASuSAv2Hybrid   || modelConfig == kMd_HFSuSAv2Hybrid ||
        modelConfig == kMd_CRPAPWSuSAv2Hybrid || modelConfig == kMd_HFPWSuSAv2Hybrid ||
        modelConfig == kMd_SuSAv2Blend)
    {
        tensor_blen = dynamic_cast<const LabFrameHadronTensorI*>
            ( fHadronTensorModel->GetTensor (tensor_pdg_crpa, tensor_type_blen) );

        if ( !tensor_blen ) 
        {
            LOG("SuSAv2QE", pWARN) << "Can't calculate final lepton polarization. Set it to zero.";
            fFinalLeptonPolarization = TVector3(0, 0, 0);
            return fFinalLeptonPolarization;
        }
    }

    if( modelConfig != kMd_SuSAv2 && modelConfig != kMd_SuSAv2Blend)
    {
        tensor_crpa = dynamic_cast<const LabFrameHadronTensorI*>
            ( fHadronTensorModel->GetTensor (tensor_pdg_crpa, tensor_type_crpa) );

        if ( !tensor_crpa ) 
        {
            LOG("SuSAv2QE", pWARN) << "Can't calculate final lepton polarization. Set it to zero.";
            fFinalLeptonPolarization = TVector3(0, 0, 0);
            return fFinalLeptonPolarization;
        }
    }

    // *****************************
    // Q_value offset calculation
    // *****************************

    // SD: The Q-Value essentially corrects q0 to account for nuclear
    // binding energy in the Valencia model but this effect is already
    // in the SuSAv2 and CRPA/HF tensors so I'll set it to 0.
    // However, if I want to scale I need to account for the altered
    // binding energy. To first order I can use the Q_value for this
    double Delta_Q_value_susa = Eb_tgt-Eb_ten_susa;
    double Delta_Q_value_crpa = Eb_tgt-Eb_ten_crpa;
    double Delta_Q_value_blen = Eb_tgt-Eb_ten_crpa;

    // Apply Qvalue relative shift if needed:
    if( fQvalueShifter ) 
    {
        // We have the option to add an additional shift on top of the binding energy correction
        // The QvalueShifter, is a relative shift to the Q_value.
        // The Q_value was already taken into account in the hadron tensor. Here we recalculate it
        // to get the right absolute shift.
        double tensor_Q_value_susa = genie::utils::mec::Qvalue(tensor_pdg_susa,probe_pdg);
        double total_Q_value_susa = tensor_Q_value_susa + Delta_Q_value_susa ;
        double Q_value_shift_susa = total_Q_value_susa * fQvalueShifter -> Shift( interaction->InitState().Tgt() ) ;

        double tensor_Q_value_crpa = genie::utils::mec::Qvalue(tensor_pdg_crpa,probe_pdg);
        double total_Q_value_crpa = tensor_Q_value_crpa + Delta_Q_value_crpa ;
        double Q_value_shift_crpa = total_Q_value_crpa * fQvalueShifter -> Shift( interaction->InitState().Tgt() ) ;

        Delta_Q_value_susa += Q_value_shift_susa;
        Delta_Q_value_crpa += Q_value_shift_crpa;
        Delta_Q_value_blen += Q_value_shift_crpa;
    }

    // Set the xsec to zero for interactions with q0,q3 outside the requested range

    if( modelConfig == kMd_SuSAv2)
    {
        double Q0min = tensor_susa->q0Min();
        double Q0max = tensor_susa->q0Max();
        double Q3min = tensor_susa->qMagMin();
        double Q3max = tensor_susa->qMagMax();
        if (Q0-Delta_Q_value_susa < Q0min || Q0-Delta_Q_value_susa > Q0max || Q3 < Q3min || Q3 > Q3max) 
        {
            LOG("SuSAv2QE", pWARN) << "Can't calculate final lepton polarization. Set it to zero.";
            fFinalLeptonPolarization = TVector3(0, 0, 0);
            return fFinalLeptonPolarization;
        }
    }
    else if ( modelConfig == kMd_CRPA   || modelConfig == kMd_HF ||
            modelConfig == kMd_CRPAPW || modelConfig == kMd_HFPW )
    {
        double Q0min = tensor_crpa->q0Min();
        double Q0max = tensor_crpa->q0Max();
        double Q3min = tensor_crpa->qMagMin();
        double Q3max = tensor_crpa->qMagMax();
        if (Q0-Delta_Q_value_crpa < Q0min || Q0-Delta_Q_value_crpa > Q0max || Q3 < Q3min || Q3 > Q3max)
        {
            LOG("SuSAv2QE", pWARN) << "Can't calculate final lepton polarization. Set it to zero.";
            fFinalLeptonPolarization = TVector3(0, 0, 0);
            return fFinalLeptonPolarization;
        }
    }
    else if ( modelConfig == kMd_SuSAv2Blend)
    {
        double Q0min = tensor_blen->q0Min();
        double Q0max = tensor_blen->q0Max();
        double Q3min = tensor_blen->qMagMin();
        double Q3max = tensor_blen->qMagMax();
        if (Q0-Delta_Q_value_blen < Q0min || Q0-Delta_Q_value_blen > Q0max || Q3 < Q3min || Q3 > Q3max) 
        {
            LOG("SuSAv2QE", pWARN) << "Can't calculate final lepton polarization. Set it to zero.";
            fFinalLeptonPolarization = TVector3(0, 0, 0);
            return fFinalLeptonPolarization;
        }
    }
    else
    { // hybrid (blending) cases. Low kinematics handled by CRPA/HF, high kinematics by blended SuSA
        double Q0min = tensor_crpa->q0Min();
        double Q0max = tensor_blen->q0Max();
        double Q3min = tensor_crpa->qMagMin();
        double Q3max = tensor_blen->qMagMax();
        if (Q0-Delta_Q_value_crpa < Q0min || Q0-Delta_Q_value_blen > Q0max || Q3 < Q3min || Q3 > Q3max)
        {
            LOG("SuSAv2QE", pWARN) << "Can't calculate final lepton polarization. Set it to zero.";
            fFinalLeptonPolarization = TVector3(0, 0, 0);
            return fFinalLeptonPolarization;
        }
    }

    double W1_susa(0), W2_susa(0), W3_susa(0), W4_susa(0), W5_susa(0)/*, W6_susa(0)*/;
    double W1_crpa(0), W2_crpa(0), W3_crpa(0), W4_crpa(0), W5_crpa(0)/*, W6_crpa(0)*/;
    double W1_blen(0), W2_blen(0), W3_blen(0), W4_blen(0), W5_blen(0)/*, W6_blen(0)*/;
    
    if ( modelConfig == kMd_SuSAv2 )
    {
        W1_susa = tensor_susa->W1(Q0 - Delta_Q_value_susa, Q3, M);
        W2_susa = tensor_susa->W2(Q0 - Delta_Q_value_susa, Q3, M);
        W3_susa = tensor_susa->W3(Q0 - Delta_Q_value_susa, Q3, M);
        W4_susa = tensor_susa->W4(Q0 - Delta_Q_value_susa, Q3, M);
        W5_susa = tensor_susa->W5(Q0 - Delta_Q_value_susa, Q3, M);
//        W6_susa = tensor_susa->W6(Q0 - Delta_Q_value_susa, Q3, M);
    }
    if ( modelConfig == kMd_CRPASuSAv2Hybrid   || modelConfig == kMd_HFSuSAv2Hybrid ||
         modelConfig == kMd_CRPAPWSuSAv2Hybrid || modelConfig == kMd_HFPWSuSAv2Hybrid ||
         modelConfig == kMd_SuSAv2Blend)
    {
        W1_blen = tensor_blen->W1(Q0 - Delta_Q_value_blen, Q3, M);
        W2_blen = tensor_blen->W2(Q0 - Delta_Q_value_blen, Q3, M);
        W3_blen = tensor_blen->W3(Q0 - Delta_Q_value_blen, Q3, M);
        W4_blen = tensor_blen->W4(Q0 - Delta_Q_value_blen, Q3, M);
        W5_blen = tensor_blen->W5(Q0 - Delta_Q_value_blen, Q3, M);
//        W6_blen = tensor_blen->W6(Q0 - Delta_Q_value_blen, Q3, M);
    }
    if( modelConfig != kMd_SuSAv2 && modelConfig != kMd_SuSAv2Blend)
    {
        W1_crpa = tensor_crpa->W1(Q0 - Delta_Q_value_crpa, Q3, M);
        W2_crpa = tensor_crpa->W2(Q0 - Delta_Q_value_crpa, Q3, M);
        W3_crpa = tensor_crpa->W3(Q0 - Delta_Q_value_crpa, Q3, M);
        W4_crpa = tensor_crpa->W4(Q0 - Delta_Q_value_crpa, Q3, M);
        W5_crpa = tensor_crpa->W5(Q0 - Delta_Q_value_crpa, Q3, M);
//        W6_crpa = tensor_crpa->W6(Q0 - Delta_Q_value_crpa, Q3, M);
    }
    
    // Apply blending if needed
    double W1(0), W2(0), W3(0), W4(0), W5(0)/*, W6(0)*/;

    if ( modelConfig == kMd_SuSAv2 )
    {      
        W1 = W1_susa;
        W2 = W2_susa;
        W3 = W3_susa;
        W4 = W4_susa;
        W5 = W5_susa;
//        W6 = W6_susa;
        
    }
    if ( modelConfig == kMd_SuSAv2Blend )
    {
        W1 = W1_blen;
        W2 = W2_blen;
        W3 = W3_blen;
        W4 = W4_blen;
        W5 = W5_blen;
//        W6 = W6_blen;
    }
    if( modelConfig == kMd_CRPA   || modelConfig == kMd_HF ||
        modelConfig == kMd_CRPAPW || modelConfig == kMd_HFPW )
    {
        W1 = W1_crpa;
        W2 = W2_crpa;
        W3 = W3_crpa;
        W4 = W4_crpa;
        W5 = W5_crpa;
//        W6 = W6_crpa;
    }
    else if( modelConfig == kMd_CRPASuSAv2Hybrid   ||
             modelConfig == kMd_HFSuSAv2Hybrid     ||
             modelConfig == kMd_CRPAPWSuSAv2Hybrid ||
             modelConfig == kMd_HFPWSuSAv2Hybrid )
    {  
        // blending cases
        if(blendMode == 1) // Linear blending in q0
            if (Q0 < q0BlendStart)
            {
                W1 = W1_crpa;
                W2 = W2_crpa;
                W3 = W3_crpa;
                W4 = W4_crpa;
                W5 = W5_crpa;
//                W6 = W6_crpa;
            }
            else if (Q0 > q0BlendEnd)
            {
                W1 = W1_blen;
                W2 = W2_blen;
                W3 = W3_blen;
                W4 = W4_blen;
                W5 = W5_blen;
//                W6 = W6_blen;
            }
            else
            {
                double SuSAFrac = (Q0 - q0BlendStart) / (q0BlendEnd - q0BlendStart);
                double CRPAFrac = 1 - SuSAFrac;
                W1 = SuSAFrac*W1_blen + CRPAFrac*W1_crpa;
                W2 = SuSAFrac*W2_blen + CRPAFrac*W2_crpa;
                W3 = SuSAFrac*W3_blen + CRPAFrac*W3_crpa;
                W4 = SuSAFrac*W4_blen + CRPAFrac*W4_crpa;
                W5 = SuSAFrac*W5_blen + CRPAFrac*W5_crpa;
//                W6 = SuSAFrac*W6_blen + CRPAFrac*W6_crpa;
            }
        else if(blendMode == 2)
        { // Exp blending in q (from Alexis)
            double phi_q = (genie::constants::kPi / 2.) * (1 - 1./(1+std::exp( (Q3 - qBlendRef)/qBlendDel)) );
            double sn2 = TMath::Sin(phi_q)*TMath::Sin(phi_q);
            double cn2 = 1 - sn2;
            W1 = sn2*W1_blen + cn2*W1_crpa;
            W2 = sn2*W2_blen + cn2*W2_crpa;
            W3 = sn2*W3_blen + cn2*W3_crpa;
            W4 = sn2*W4_blen + cn2*W4_crpa;
            W5 = sn2*W5_blen + cn2*W5_crpa;
//            W6 = sn2*W6_blen + cn2*W6_crpa;
        }
    }

    double p[4], q[4], epq[4][4], k[4], l[4], s[4], eskl[4];
    std::complex<double> jp[4], jm[4];
    
    p[0] = inNucleonMomOnShell.E();
    p[1] = inNucleonMomOnShell.Px();
    p[2] = inNucleonMomOnShell.Py();
    p[3] = inNucleonMomOnShell.Pz();
  
    q[0] = qTildeP4.E();
    q[1] = qTildeP4.Px();
    q[2] = qTildeP4.Py();
    q[3] = qTildeP4.Pz();
  
    k[0] = neutrinoMom.E();
    k[1] = -neutrinoMom.Px();
    k[2] = -neutrinoMom.Py();
    k[3] = -neutrinoMom.Pz();
  
    l[0] = leptonMom.E();
    l[1] = -leptonMom.Px();
    l[2] = -leptonMom.Py();
    l[3] = -leptonMom.Pz();
  
    s[0] = leptonMom.P()/ml;
    s[1] = -leptonMom.Vect().Unit().X()*leptonMom.E()/ml;
    s[2] = -leptonMom.Vect().Unit().Y()*leptonMom.E()/ml;
    s[3] = -leptonMom.Vect().Unit().Z()*leptonMom.E()/ml;
  
    // epsilon^\alpha\beta\gamma\delta p_\gamma q_\delta
    for (int a = 0; a < 4; a++)
    {
        for (int b = 0; b < 4; b++)
        {
            epq[a][b] = 0;
            if (b == a) continue;
            for (int g = 0; g < 4; g++)
            {
                if (g == b || g == a) continue;
                for (int d = 0; d < 4; d++)
                {
                    if (d == g || d == b || d == a) continue;
                    epq[a][b] += e(a,b,g,d)*(a == 0?1:-1)*(b == 0?1:-1)*p[g]*q[d];
                }
            }
        }
    }
    
    // epsilon_\alpha\beta\gamma\delta s^\beta k^\gamma l^\delta
    for (int a = 0; a < 4; a++)
    {
        eskl[a] = 0;
        for (int b = 0; b < 4; b++)
        {
            if (b == a) continue;
            for (int g = 0; g < 4; g++)
            {
                if (g == b || g == a) continue;
                for (int d = 0; d < 4; d++)
                {
                    if (d == g || d == b || d == a) continue;
                    double sb = s[b]*(b == 0?1:-1);
                    double kg = k[g]*(g == 0?1:-1);
                    double ld = l[d]*(d == 0?1:-1);
                    eskl[a] += e(a,b,g,d)*sb*kg*ld;
                }
            }
        }
    }
        
    double kl = k[0]*l[0] - k[1]*l[1] - k[2]*l[2] - k[3]*l[3];
    double ks = k[0]*s[0] - k[1]*s[1] - k[2]*s[2] - k[3]*s[3];
            
    for (int a = 0; a < 4; a++)
    {
        if (is_neutrino)
        {
            jp[a] =  (l[a]*ks - s[a]*kl - 1i*eskl[a] + ml*k[a])/sqrt(kl + ml*ks);   //jp_\alpha
            jm[a] = (-l[a]*ks + s[a]*kl + 1i*eskl[a] + ml*k[a])/sqrt(kl - ml*ks);   //jm_\alpha
        }
        else
        {
            jp[a] =  (l[a]*ks - s[a]*kl + 1i*eskl[a] - ml*k[a])/sqrt(kl - ml*ks);   //jp_\alpha
            jm[a] =  (l[a]*ks - s[a]*kl + 1i*eskl[a] + ml*k[a])/sqrt(kl + ml*ks);   //jm_\alpha
        }
    }
    
    
    //Additional constants and variables
    std::complex<double> Wmunu, Wnumu, LWpp(0, 0), LWpm(0, 0), LWmp(0, 0), LWmm(0, 0);
    for(int mu = 0; mu < 4; mu++)
    {
        for(int nu = mu;nu < 4; nu++)
        {
            double Wreal = -g(mu,nu)*W1 + p[mu]*p[nu]*W2/M2 + q[mu]*q[nu]*W4/M2 + (p[mu]*q[nu] + q[mu]*p[nu])*W5/2/M2;
            double Wimag = epq[mu][nu]*W3/2/M2;
            Wmunu = Wreal - 1i*Wimag;  // W^\mu\nu
            LWpp += jp[mu]*std::conj(jp[nu])*Wmunu; // Lpp_\mu\nu*W^\mu\nu
            LWpm += jp[mu]*std::conj(jm[nu])*Wmunu; // Lpm_\mu\nu*W^\mu\nu
            LWmp += jm[mu]*std::conj(jp[nu])*Wmunu; // Lmp_\mu\nu*W^\mu\nu
            LWmm += jm[mu]*std::conj(jm[nu])*Wmunu; // Lmm_\mu\nu*W^\mu\nu
            if (mu != nu)
            {
                Wnumu = Wreal + 1i*Wimag;
                LWpp += jp[nu]*std::conj(jp[mu])*Wnumu; // Lpp_\mu\nu*W^\mu\nu
                LWpm += jp[nu]*std::conj(jm[mu])*Wnumu; // Lpm_\mu\nu*W^\mu\nu
                LWmp += jm[nu]*std::conj(jp[mu])*Wnumu; // Lmp_\mu\nu*W^\mu\nu
                LWmm += jm[nu]*std::conj(jm[mu])*Wnumu; // Lmm_\mu\nu*W^\mu\nu
            }
        }
    }
    std::complex<double> LWppmm = LWpp + LWmm;
    std::complex<double> rhopp = LWpp/LWppmm;
    std::complex<double> rhopm = LWpm/LWppmm;
    std::complex<double> rhomp = LWmp/LWppmm;
    std::complex<double> rhomm = LWmm/LWppmm;
    double PL = std::real(rhopp - rhomm);
    double PP = std::real(rhopm + rhomp);
    double PT = std::imag(rhomp - rhopm);
    
    TVector3 neutrinoMom3 = neutrinoMom.Vect();                                          
    TVector3 leptonMom3 = leptonMom.Vect();
    TVector3 Pz = leptonMom3.Unit();
    TVector3 Px = neutrinoMom3.Cross(leptonMom3).Unit();
    TVector3 Py = Pz.Cross(Px);
    TVector3 pol = PT*Px + PP*Py + PL*Pz;
    fFinalLeptonPolarization = pol;
    
    std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    std::cout << "PL = " << PL << ", PT = " << PT << ", PP = " << PP << "\n";
    std::cout << fFinalLeptonPolarization.Mag() << "\n";
    std::cout << "SU@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << std::endl;
    
    return fFinalLeptonPolarization;
}
//____________________________________________________________________________
inline int SuSAv2QELPXSec::g(int a, int b) const
{
    return (a==b)*(2*(a==0) - 1);
}
//____________________________________________________________________________
inline int SuSAv2QELPXSec::e(int a, int b, int c, int d) const
{
    return (b - a)*(c - a)*(d - a)*(c - b)*(d - b)*(d - c)/12;
}
//____________________________________________________________________________
