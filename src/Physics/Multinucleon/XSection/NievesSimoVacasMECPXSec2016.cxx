//_________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Original code contributed by J. Schwehr, D. Cherdack, R. Gran
 Substantial code refactorizations by the core GENIE group.
*/
//_________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/HadronTensors/NievesMECHadronTensorModel.h"
#include "Physics/HadronTensors/ValenciaHadronTensorI.h"
#include "Physics/Multinucleon/XSection/NievesSimoVacasMECPXSec2016.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
// For q3q0 scan verification and RW histos:
#include <TH1F.h>
#include <TH3D.h>
#include <TFile.h>
#include "Framework/Conventions/Constants.h"

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

  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>(
    fHadronTensorModel->GetTensor(tensor_pdg, genie::kHT_MEC_FullAll) );

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

    const ValenciaHadronTensorI* tensor_delta_all
      = dynamic_cast<const ValenciaHadronTensorI*>(
      fHadronTensorModel->GetTensor(tensor_pdg, genie::kHT_MEC_DeltaAll) );

    if ( !tensor_delta_all ) {
      LOG("NievesSimoVacasMEC", pWARN) << "Failed to load a \"DeltaAll\""
        << " hadronic tensor for nuclide " << tensor_pdg;
      return 0.;
    }

    const ValenciaHadronTensorI* tensor_delta_pn
      = dynamic_cast<const ValenciaHadronTensorI*>(
      fHadronTensorModel->GetTensor(tensor_pdg, genie::kHT_MEC_Deltapn) );

    if ( !tensor_delta_pn ) {
      LOG("NievesSimoVacasMEC", pWARN) << "Failed to load a \"Deltapn\""
        << " hadronic tensor for nuclide " << tensor_pdg;
      return 0.;
    }

    xsec_all = tensor_delta_all->dSigma_dT_dCosTheta(interaction, Q_value);
    xsec_pn = tensor_delta_pn->dSigma_dT_dCosTheta(interaction, Q_value);

  }
  else {

    const ValenciaHadronTensorI* tensor_full_all
      = dynamic_cast<const ValenciaHadronTensorI*>(
      fHadronTensorModel->GetTensor(tensor_pdg, genie::kHT_MEC_FullAll) );

    if ( !tensor_full_all ) {
      LOG("NievesSimoVacasMEC", pWARN) << "Failed to load a \"FullAll\""
        << " hadronic tensor for nuclide " << tensor_pdg;
      return 0.;
    }

    const ValenciaHadronTensorI* tensor_full_pn
      = dynamic_cast<const ValenciaHadronTensorI*>(
      fHadronTensorModel->GetTensor(tensor_pdg, genie::kHT_MEC_Fullpn) );

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

	fHadronTensorModel = dynamic_cast<const HadronTensorModelI *> (
          this->SubAlg("HadronTensorAlg") );
        assert( fHadronTensorModel );

	fXSecIntegrator =
        dynamic_cast<const XSecIntegratorI *> (
          this->SubAlg("NumericalIntegrationAlg"));
        assert(fXSecIntegrator);

  //uncomment this to make histo output on a 2D grid for validation
  //this->Scanq0q3();
  //this->Scanq0q3_np();

  //uncomment this to make 3D histograms to be used for RW-ing
  //this->buildRWhistos();
}
//_________________________________________________________________________
void NievesSimoVacasMECPXSec2016::Scanq0q3(void)
{
  std::cout << "Calling NievesSimoVacasMECPXSec2016::Scanq0q3" << std::endl;

  // Runs after loading the tensors.  But all this is hard coded, so you need a private build.
  // Can gmkspl | grep dsdq0q3 > output.dat

  int nu_pdg = 14;
  int target_pdg = kPdgTgtC12;
  double Enu = 3.0;
  double m_probe = 0.0;
  double Ml = PDGLibrary::Instance()->Find(13)->Mass();;
  //double Ml = kElectronMass;
  double Tl = 0.0;
  double costhl = 0.0;
  int tensor_pdg = kPdgTgtC12;
  HadronTensorType_t tensor_type = kHT_MEC_FullAll;

  int nbins = 500;
  double upper = 2.0;
  double step = upper/(double)nbins;
  double addhalfstep = step/2.0;

  double corr_fact=0.0;
  double myq = 0.0;
  double myw = 0.0;

  //double Q_value = 0.0;
  double Q_value = genie::utils::mec::Qvalue(target_pdg, nu_pdg);

  TFile* nievesq0q3 = new TFile("nievesq0q3.root","RECREATE");
  TH2D*  nievesq0q3Hist = new TH2D("nievesq0q3", "nievesq0q3", nbins, 0, upper, nbins, 0, upper);
  TH2D*  nievesq0q3Hist_GU = new TH2D("nievesq0q3_GU", "nievesq0q3_GU", nbins, 0, upper, nbins, 0, upper);
  TH2D*  nievesCthTlHist = new TH2D("nievesCthTlHist", "nievesCthTlHist", 100, -1, 1, 100, 0, upper);

  TH2D*  nievesW00Hist = new TH2D("nievesW00Hist", "nievesW00Hist", nbins, 0, upper, nbins, 0, upper); //tt
  TH2D*  nievesW03Hist = new TH2D("nievesW03Hist", "nievesW03Hist", nbins, 0, upper, nbins, 0, upper); //tz
  TH2D*  nievesW11Hist = new TH2D("nievesW11Hist", "nievesW11Hist", nbins, 0, upper, nbins, 0, upper); //xx
  TH2D*  nievesW12Hist = new TH2D("nievesW12Hist", "nievesW12Hist", nbins, 0, upper, nbins, 0, upper); //xy
  TH2D*  nievesW33Hist = new TH2D("nievesW33Hist", "nievesW33Hist", nbins, 0, upper, nbins, 0, upper); //zz


  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>(
    fHadronTensorModel->GetTensor(target_pdg, tensor_type) );


  double Q0min = tensor->q0Min();
  double Q0max = tensor->q0Max();
  double Q3min = tensor->qMagMin();
  double Q3max = tensor->qMagMax();

  std::cout << "Scanq0q3: q0 min / max is " << Q0min << " / " << Q0max << std::endl;
  std::cout << "Scanq0q3: q3 min / max is " << Q3min << " / " << Q3max << std::endl;

  for(int iq = 0; iq<nbins; iq++){
    for(int iw = 0; iw<nbins; iw++){
      myq = addhalfstep + step * (double)iq;
      myw = addhalfstep + step * (double)iw;
      int myBin_qw = nievesq0q3Hist->FindBin(myq,myw);


      genie::utils::mec::GetTlCostlFromq0q3(myw, myq, Enu, Ml, Tl, costhl);

      double xsec = 0.0;
      if(myq > Q3max || myq < Q3min || myw > Q0max || myw < Q0min){
        xsec = 0.0;
      }
      else if (Tl < 0.0){
        xsec = 0.0;
      }
      else {
        xsec = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);

        double w00 = (tensor->tt(myw,myq)).real();
        double w03 = (tensor->tz(myw,myq)).real();
        double w11 = (tensor->xx(myw,myq)).real();
        double w12 = -1.0*(tensor->xy(myw,myq)).imag();
        double w33 = (tensor->zz(myw,myq)).real();

        nievesW00Hist->SetBinContent(myBin_qw, w00);
        nievesW03Hist->SetBinContent(myBin_qw, w03);
        nievesW11Hist->SetBinContent(myBin_qw, w11);
        nievesW12Hist->SetBinContent(myBin_qw, w12);
        nievesW33Hist->SetBinContent(myBin_qw, w33);

      }

      nievesq0q3Hist_GU->SetBinContent(myBin_qw, xsec);

      // Apply the squared CKM matrix element Vud as a correction factor
      // (not included in the tabulated tensor element values)
      //xsec *= fVud2;

      // Apply given overall scaling factor
      //xsec *= fXSecScale;

      xsec = xsec / (1.0E-38 * units::cm2);  // output in something other than genie units.

      //SD 04/12/17 - apply SuSAv2 correction factor (cos^2theta_c):
      //corr_fact=0.9471182*1000.;

      nievesq0q3Hist->SetBinContent(myBin_qw, xsec);

      int myBin_CthTl = nievesCthTlHist->FindBin(costhl,Tl);
      nievesCthTlHist->SetBinContent(myBin_CthTl, xsec);

    }
  }
  nievesq0q3Hist->Write();
  nievesq0q3Hist_GU->Write();
  nievesCthTlHist->Write();
  nievesW00Hist->Write();
  nievesW03Hist->Write();
  nievesW11Hist->Write();
  nievesW12Hist->Write();
  nievesW33Hist->Write();
  nievesq0q3Hist->SaveAs("nievesq0q3.png");
  nievesq0q3->Close();
}
//_________________________________________________________________________
void NievesSimoVacasMECPXSec2016::Scanq0q3_np(void)
{
  std::cout << "Calling NievesSimoVacasMECPXSec2016::Scanq0q3_np" << std::endl;

  // Runs after loading the tensors.  But all this is hard coded, so you need a private build.
  // Can gmkspl | grep dsdq0q3 > output.dat

  int nu_pdg = 14;
  int target_pdg = kPdgTgtC12;
  double Enu = 3.0;
  double m_probe = 0.0;
  double Ml = PDGLibrary::Instance()->Find(13)->Mass();;
  //double Ml = kElectronMass;
  double Tl = 0.0;
  double costhl = 0.0;
  int tensor_pdg = kPdgTgtC12;
  HadronTensorType_t tensor_type = kHT_MEC_Fullpn;

  int nbins = 500;
  double upper = 2.0;
  double step = upper/(double)nbins;
  double addhalfstep = step/2.0;

  double corr_fact=0.0;
  double myq = 0.0;
  double myw = 0.0;

  //double Q_value = 0.0;
  double Q_value = genie::utils::mec::Qvalue(target_pdg, nu_pdg);

  TFile* nievesq0q3 = new TFile("nievesq0q3_np.root","RECREATE");
  TH2D*  nievesq0q3Hist = new TH2D("nievesq0q3", "nievesq0q3", nbins, 0, upper, nbins, 0, upper);
  TH2D*  nievesq0q3Hist_GU = new TH2D("nievesq0q3_GU", "nievesq0q3_GU", nbins, 0, upper, nbins, 0, upper);
  TH2D*  nievesCthTlHist = new TH2D("nievesCthTlHist", "nievesCthTlHist", 100, -1, 1, 100, 0, upper);

  TH2D*  nievesW00Hist = new TH2D("nievesW00Hist", "nievesW00Hist", nbins, 0, upper, nbins, 0, upper); //tt
  TH2D*  nievesW03Hist = new TH2D("nievesW03Hist", "nievesW03Hist", nbins, 0, upper, nbins, 0, upper); //tz
  TH2D*  nievesW11Hist = new TH2D("nievesW11Hist", "nievesW11Hist", nbins, 0, upper, nbins, 0, upper); //xx
  TH2D*  nievesW12Hist = new TH2D("nievesW12Hist", "nievesW12Hist", nbins, 0, upper, nbins, 0, upper); //xy
  TH2D*  nievesW33Hist = new TH2D("nievesW33Hist", "nievesW33Hist", nbins, 0, upper, nbins, 0, upper); //zz


  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>(
    fHadronTensorModel->GetTensor(target_pdg, tensor_type) );

  double Q0min = tensor->q0Min();
  double Q0max = tensor->q0Max();
  double Q3min = tensor->qMagMin();
  double Q3max = tensor->qMagMax();

  std::cout << "Scanq0q3: q0 min / max is " << Q0min << " / " << Q0max << std::endl;
  std::cout << "Scanq0q3: q3 min / max is " << Q3min << " / " << Q3max << std::endl;

  for(int iq = 0; iq<nbins; iq++){
    for(int iw = 0; iw<nbins; iw++){
      myq = addhalfstep + step * (double)iq;
      myw = addhalfstep + step * (double)iw;
      int myBin_qw = nievesq0q3Hist->FindBin(myq,myw);


      genie::utils::mec::GetTlCostlFromq0q3(myw, myq, Enu, Ml, Tl, costhl);

      double xsec = 0.0;
      if(myq > Q3max || myq < Q3min || myw > Q0max || myw < Q0min){
        xsec = 0.0;
      }
      else if (Tl < 0.0){
        xsec = 0.0;
      }
      else {
        xsec = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);

        double w00 = (tensor->tt(myw,myq)).real();
        double w03 = (tensor->tz(myw,myq)).real();
        double w11 = (tensor->xx(myw,myq)).real();
        double w12 = -1.0*(tensor->xy(myw,myq)).imag();
        double w33 = (tensor->zz(myw,myq)).real();

        nievesW00Hist->SetBinContent(myBin_qw, w00);
        nievesW03Hist->SetBinContent(myBin_qw, w03);
        nievesW11Hist->SetBinContent(myBin_qw, w11);
        nievesW12Hist->SetBinContent(myBin_qw, w12);
        nievesW33Hist->SetBinContent(myBin_qw, w33);

      }

      nievesq0q3Hist_GU->SetBinContent(myBin_qw, xsec);

      // Apply the squared CKM matrix element Vud as a correction factor
      // (not included in the tabulated tensor element values)
      //xsec *= fVud2;

      // Apply given overall scaling factor
      //xsec *= fXSecScale;

      xsec = xsec / (1.0E-38 * units::cm2);  // output in something other than genie units.

      //SD 04/12/17 - apply SuSAv2 correction factor (cos^2theta_c):
      //corr_fact=0.9471182*1000.;

      nievesq0q3Hist->SetBinContent(myBin_qw, xsec);

      int myBin_CthTl = nievesCthTlHist->FindBin(costhl,Tl);
      nievesCthTlHist->SetBinContent(myBin_CthTl, xsec);

    }
  }
  nievesq0q3Hist->Write();
  nievesq0q3Hist_GU->Write();
  nievesCthTlHist->Write();
  nievesW00Hist->Write();
  nievesW03Hist->Write();
  nievesW11Hist->Write();
  nievesW12Hist->Write();
  nievesW33Hist->Write();
  nievesq0q3Hist->SaveAs("nievesq0q3.png");
  nievesq0q3->Close();
}

//_________________________________
void NievesSimoVacasMECPXSec2016::buildRWhistos(void)
{
  std::cout << "Calling NievesSimoVacasMECPXSec2016::buildRWhistos" << std::endl;

  // Runs after loading the tensors.  But all this is hard coded, so you need a private build.

  int nu_pdg = 14;
  int target_pdg = kPdgTgtC12;
  double Enu = 3.0;
  double m_probe = 0.0;
  double Ml = PDGLibrary::Instance()->Find(13)->Mass();;
  //double Ml = kElectronMass;
  double Tl = 0.0;
  double pmu = 0.0;
  double costhl = 0.0;
  int tensor_pdg = kPdgTgtC12;
  HadronTensorType_t tensor_type = kHT_MEC_FullAll;

  int nbins = 100;
  double upper = 2.0;
  double step = upper/(double)nbins;
  double addhalfstep = step/2.0;

  double upperEnu = 10.0;
  double stepEnu = upperEnu/(double)nbins;
  double addhalfstepEnu = stepEnu/2.0;

  double corr_fact=0.0;
  double myq = 0.0;
  double myw = 0.0;
  double myenu = 0.0;

  double Q_value = 0.0;

  TFile* nievesq0q3 = new TFile("nievesMECRWHists.root","RECREATE");
  TH3D* nievesEnuq0q3Hist = new TH3D("nievesElq0q3", "nievesElq0q3", nbins, 0, upperEnu, nbins, 0, upper, nbins, 0, upper);
  TH3D* nievesEnuCthTlHist = new TH3D("nievesElCthTlHist", "nievesElCthTlHist", nbins, 0, upperEnu, 100, -1, 1, 100, 0, upper);



  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>(
    fHadronTensorModel->GetTensor(target_pdg, tensor_type) );

  if(!tensor){
    std::cout << "ERROR!!!! Tensor is NULL" << std::endl;
  }

  double Q0min = tensor->q0Min();
  double Q0max = tensor->q0Max();
  double Q3min = tensor->qMagMin();
  double Q3max = tensor->qMagMax();

  std::cout << "buildRWhistos: q0 min / max is " << Q0min << " / " << Q0max << std::endl;
  std::cout << "buildRWhistos: q3 min / max is " << Q3min << " / " << Q3max << std::endl;


  for(int ienu = 0; ienu<nbins; ienu++){
    for(int iq = 0; iq<nbins; iq++){
      for(int iw = 0; iw<nbins; iw++){
        myenu = addhalfstepEnu + stepEnu * (double)ienu;
        myq = addhalfstep + step * (double)iq;
        myw = addhalfstep + step * (double)iw;
        int myBin_qw = nievesEnuq0q3Hist->FindBin(myenu,myq,myw);

        if(ienu%50 && iq==0 && iw==0) std::cout << "buildRWhistos: energy is " << myenu <<  ", max is " << upperEnu << std::endl;

        genie::utils::mec::GetTlCostlFromq0q3(myw, myq, Enu, Ml, Tl, costhl);

        double xsec = 0.0;
        if(myq > Q3max || myq < Q3min || myw > Q0max || myw < Q0min){
          xsec = 0.0;
        }
        else if (Tl < 0.0){
          xsec = 0.0;
        }
        else {
          xsec = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
        }

        // Apply given overall scaling factor
        xsec *= fXSecScale;

        // output in something other than genie units.
        xsec = xsec / (1.0E-38 * units::cm2);

        //Find bin and fill:
        nievesEnuq0q3Hist->SetBinContent(myBin_qw, xsec);

        int myBin_CthTl = nievesEnuCthTlHist->FindBin(myenu, costhl,Tl);
        nievesEnuCthTlHist->SetBinContent(myBin_CthTl, xsec);
        //nievesq0q3Hist->Fill(myq,myw,xsec);

      }
    }
  }

  nievesEnuq0q3Hist->Write();
  nievesEnuCthTlHist->Write();

  //nievesq0q3Hist->SaveAs("nievesq0q3.png");
  nievesq0q3->Close();
}
