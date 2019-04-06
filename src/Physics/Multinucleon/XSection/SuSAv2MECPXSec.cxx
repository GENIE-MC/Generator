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
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"

// For q3q0 scan verification and RW histos:
#include <TH1F.h>
#include <TH3D.h>
#include <TFile.h>
#include "Framework/Conventions/Constants.h"

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
  int tensor_pdg = target_pdg;
  int A_request = pdg::IonPdgCodeToA(target_pdg);
  int Z_request = pdg::IonPdgCodeToZ(target_pdg);
  bool need_to_scale = false;

  HadronTensorType_t tensor_type = kHT_Undefined;
  HadronTensorType_t pn_tensor_type = kHT_Undefined;
  if ( pdg::IsNeutrino(probe_pdg) || pdg::IsAntiNeutrino(probe_pdg) ) {
    tensor_type = kHT_MEC_FullAll;
    pn_tensor_type = kHT_MEC_Fullpn; 
    //tensor_type = kHT_MEC_FullAll_wImag;
    //pn_tensor_type = kHT_MEC_FullAll_wImag; 
  }
  else {
    // If the probe is not a neutrino, assume that it's an electron
    // For the moment all electron interactions are pp final state
    tensor_type = kHT_MEC_EM;
    pn_tensor_type = kHT_MEC_EM;
  }

  /// \todo Add more hadron tensors so this scaling is not so terrible
  if ( A_request == 4 && Z_request == 2 ) {
    tensor_pdg = kPdgTgtC12;
    // This is for helium 4, but use carbon tensor
    // the use of nuclear density parameterization is suspicious
    // but some users (MINERvA) need something not nothing.
  }
  else if (A_request < 9) {
    // refuse to do D, T, He3, Li, and some Be, B
    // actually it would work technically, maybe except D, T
    MAXLOG("SuSAv2MEC", pWARN, 10) << "Asked to scale to D through B, this will not work";
    return 0;
  }
  else if (A_request >= 9 && A_request < 15) {
    tensor_pdg = kPdgTgtC12;
  }
  else if(A_request >= 15 && A_request < 22) {
      tensor_pdg = kPdgTgtO16;    
  }
  else if(A_request >= 22) {
      tensor_pdg = kPdgTgtO16;
      LOG("SuSAv2MEC", pWARN) << "Dangerous scaling from " << target_pdg << " to " << tensor_pdg;    
      LOG("SuSAv2MEC", pWARN) << "Hadron tensors for such heavy elements not ready yet";    
      LOG("SuSAv2MEC", pWARN) << "Proceed at your own peril!";    
  }


  if(tensor_pdg != target_pdg) need_to_scale = true;

  // The SuSAv2-MEC hadron tensors are defined using the same conventions
  // as the Valencia MEC model, so we can use the same sort of tensor
  // object to describe them.
  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>( htp.GetTensor(tensor_pdg,
    tensor_type, fHadronTensorTableName) );

  const ValenciaHadronTensorI* tensor_pn
    = dynamic_cast<const ValenciaHadronTensorI*>( htp.GetTensor(tensor_pdg,
    pn_tensor_type, fHadronTensorTableName) );


  /// \todo add the different pair configurations for e-scattering

  // If retrieving the tensor failed, complain and return zero
  if ( !tensor ) {
    LOG("SuSAv2MEC", pWARN) << "Failed to load a hadronic tensor for the"
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
  // SD: The Q-Value essentially corrects q0 to account for nuclear
  // binding energy in the Valencia model but I this effect is already
  // in Guille's tensors so I'll set it to 0.

  int nu_pdg = interaction->InitState().ProbePdg();
  //double Q_value = genie::utils::mec::Qvalue(target_pdg, nu_pdg);
  double Q_value = 0.;

  // By default, we will compute the full cross-section. If a {p,n} hit 
  // dinucleon was set we will calculate the cross-section for that 
  // component only

  bool pn = (interaction->InitState().Tgt().HitNucPdg() == kPdgClusterNP);

  // Compute the cross section using the hadron tensor
  double xsec_all = tensor->dSigma_dT_dCosTheta_rosenbluth(interaction, Q_value);
  double xsec_pn = tensor_pn->dSigma_dT_dCosTheta_rosenbluth(interaction, Q_value);
  //double xsec_all = tensor->dSigma_dT_dCosTheta(interaction, Q_value);
  //double xsec_pn = tensor_pn->dSigma_dT_dCosTheta(interaction, Q_value);

  //hadron tensor precision can sometimes lead to 0 xsec_pn but finite xsec
  //seems to cause issues downstream ...
  if(xsec_pn==0) xsec_pn = 0.00001*xsec_all;

  // Choose the right kind of cross section ("all" or "pn") to return
  // based on whether a {p, n} dinucleon was hit
  double xsec = (pn) ? xsec_pn : xsec_all;

  
  // This scaling should be okay-ish for the total xsec, but it misses 
  // the energy shift. To get this we should really just build releveant 
  // hadron tensors but there may be some ways to approximate it.
  // For more details see Guille's thesis: https://idus.us.es/xmlui/handle/11441/74826
  if ( need_to_scale ) {
    FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
    const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
    double KF_tgt = kft->FindClosestKF(target_pdg, kPdgProton);
    double KF_ten = kft->FindClosestKF(tensor_pdg, kPdgProton);
    LOG("SuSAv2MEC", pDEBUG) << "KF_tgt = " << KF_tgt;
    LOG("SuSAv2MEC", pDEBUG) << "KF_ten = " << KF_ten;
    double A_ten  = pdg::IonPdgCodeToA(tensor_pdg);
    double scaleFact = (A_ten/A_request)*(KF_ten/KF_tgt)*(KF_ten/KF_tgt);
    xsec *= scaleFact;
  }

  // Apply the squared CKM matrix element Vud as a correction factor
  // (not included in the tabulated tensor element values)
  // NOT NEEDED IF USING THE ROSENBLUTH FORMALISM
  //xsec *= fVud2;

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
  // No longer needed as using Rosenbluth formalism
  //double Vud;
  //GetParam("CKM-Vud", Vud);
  //fVud2 = std::pow(Vud, 2);

  // Cross section scaling factor
  GetParamDef("MEC-XSecScale", fXSecScale, 1.) ;

  // Name of the hadron tensor table to use for this model
  GetParamDef("HadronTensorTableName", fHadronTensorTableName,
    std::string("SuSAv2_MEC"));

  fXSecIntegrator = dynamic_cast<const XSecIntegratorI*> (
    this->SubAlg("NumericalIntegrationAlg"));
  assert(fXSecIntegrator);

  //Fermi momentum tables for scaling
  this->GetParam( "FermiMomentumTable", fKFTable);

  //uncomment this to make histo output on a 2D grid for validation
  //this->Scanq0q3();
  //this->Scanq0q3_np();
  //this->Scanq0q3_electron();

  //uncomment this to make 3D histograms to be used for RW-ing
  //this->buildRWhistos();
}
//_________________________________________________________________________
void SuSAv2MECPXSec::Scanq0q3(void)
{
  std::cout << "Calling SuSAv2MECPXSec::Scanq0q3" << std::endl;


  // Get the hadron tensor pool
  HadronTensorPool& htp = HadronTensorPool::Instance();

  // Runs after loading the tensors.  But all this is hard coded, so you need a private build.
  // Can gmkspl | grep dsdq0q3 > output.dat
  
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
  HadronTensorType_t tensor_type = kHT_MEC_FullAll; // kHT_MEC_FullAll_wImag;

  int nbins = 500;
  double upper = 2.0;
  double step = upper/(double)nbins;
  double addhalfstep = step/2.0;

  double corr_fact=0.0;
  double myq = 0.0;
  double myw = 0.0;

  double Q_value = 0.0;

  TFile* susaq0q3 = new TFile("susaq0q3.root","RECREATE");
  TH2D* susaq0q3Hist = new TH2D("SuSAv2q0q3", "SuSAv2q0q3", nbins, 0, upper, nbins, 0, upper);
  TH2D* susaCthTlHist = new TH2D("susaCthTlHist", "susaCthTlHist", 100, -1, 1, 100, 0, upper);

  TH2D*  susaW00Hist = new TH2D("susaW00Hist", "susaW00Hist", nbins, 0, upper, nbins, 0, upper); //tt
  TH2D*  susaW03Hist = new TH2D("susaW03Hist", "susaW03Hist", nbins, 0, upper, nbins, 0, upper); //tz
  TH2D*  susaW11Hist = new TH2D("susaW11Hist", "susaW11Hist", nbins, 0, upper, nbins, 0, upper); //xx
  TH2D*  susaW12Hist = new TH2D("susaW12Hist", "susaW12Hist", nbins, 0, upper, nbins, 0, upper); //xy
  TH2D*  susaW33Hist = new TH2D("susaW33Hist", "susaW33Hist", nbins, 0, upper, nbins, 0, upper); //zz

  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>( htp.GetTensor(target_pdg,
    tensor_type, "SuSAv2_MEC") );

  if(!tensor){
    std::cout << "ERROR!!!! Tensor is NULL" << std::endl;
  }

  double Q0min = tensor->q0Min();
  double Q0max = tensor->q0Max();
  double Q3min = tensor->qMagMin();
  double Q3max = tensor->qMagMax();

  std::cout << "Scanq0q3: q0 min / max is " << Q0min << " / " << Q0max << std::endl;
  std::cout << "Scanq0q3: q3 min / max is " << Q3min << " / " << Q3max << std::endl;

  std::cout << "W00 for q0=q3=100MeV  is: " << tensor->tt(0.1, 0.1) << std::endl;
  std::cout << "W00 for q0=q3=300MeV  is: " << tensor->tt(0.3, 0.3) << std::endl;
  std::cout << "W00 for q0=q3=500MeV  is: " << tensor->tt(0.5, 0.5) << std::endl;
  std::cout << "W00 for q0=q3=700MeV  is: " << tensor->tt(0.7, 0.7) << std::endl;
  std::cout << "W00 for q0=q3=900MeV  is: " << tensor->tt(0.9, 0.9) << std::endl;
  std::cout << "W00 for q0=q3=1100MeV is: " << tensor->tt(1.1, 1.1) << std::endl;
  std::cout << "W00 for q0=q3=1300MeV is: " << tensor->tt(1.3, 1.3) << std::endl << std::endl;

  std::cout << "W11 for q0=q3=100MeV  is: " << tensor->xx(0.1, 0.1) << std::endl;
  std::cout << "W11 for q0=q3=300MeV  is: " << tensor->xx(0.3, 0.3) << std::endl;
  std::cout << "W11 for q0=q3=500MeV  is: " << tensor->xx(0.5, 0.5) << std::endl;
  std::cout << "W11 for q0=q3=700MeV  is: " << tensor->xx(0.7, 0.7) << std::endl;
  std::cout << "W11 for q0=q3=900MeV  is: " << tensor->xx(0.9, 0.9) << std::endl;
  std::cout << "W11 for q0=q3=1100MeV is: " << tensor->xx(1.1, 1.1) << std::endl;
  std::cout << "W11 for q0=q3=1300MeV is: " << tensor->xx(1.3, 1.3) << std::endl << std::endl;

  std::cout << "W11 for q0=500MeV, q3=550MeV  is: " << tensor->xx(0.5, 0.55) << std::endl;
  std::cout << "W11 for q0=550MeV, q3=500MeV  is: " << tensor->xx(0.55, 0.5) << std::endl << std::endl;


  std::cout << "W11 for q0=300MeV, q3=400MeV  is: " << tensor->xx(0.3, 0.4) << std::endl;
  std::cout << "W00 for q0=300MeV, q3=400MeV  is: " << tensor->tt(0.3, 0.4) << std::endl;

  std::cout << "W11 for q0=300MeV, q3=500MeV  is: " << tensor->xx(0.3, 0.5) << std::endl;
  std::cout << "W00 for q0=300MeV, q3=500MeV  is: " << tensor->tt(0.3, 0.5) << std::endl;

  std::cout << "W11 for q0=350MeV, q3=400MeV  is: " << tensor->xx(0.35, 0.4) << std::endl;
  std::cout << "W00 for q0=350MeV, q3=400MeV  is: " << tensor->tt(0.35, 0.4) << std::endl;

  // Fixed point debugging
  double xsec_t = 0;
  double xsec_r = 0;
  double cont_t = 0;
  double w00 = 0;
  double w03 = 0;
  double w11 = 0;
  double w12 = 0;
  double w33 = 0;
  std::cout << "Calculate some fixed point xsec: " << std::endl;
  std::cout << "Enu, pmu, Tl, costhl, w, q: xsec, xsec_r, contraction" << std::endl;
  std::cout << "       HTele: w00, w03, w11, w12, w33" << std::endl;
  Enu = 1.0;
  pmu = 0.6;
  costhl = 0.8;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  //cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;


  Enu = 0.75;
  costhl = 0.5;
  pmu = 0.6;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.75;
  costhl = 0.5;
  pmu = 0.35;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 3.0;
  costhl = 0.95;
  pmu = 2.4;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 1.5;
  costhl = 0.9;
  pmu = 1.05;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.65;
  costhl = 0.902;
  pmu = 0.21065;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;


  Enu = 0.65;
  costhl = 0.902;
  pmu = 0.4;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.6;
  costhl = 0.99;
  pmu = 0.18;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.6;
  costhl = 0.99;
  pmu = 0.35;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.6;
  costhl = 0.99;
  pmu = 0.28;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 3.0; //final Enu to use


  for(int iq = 0; iq<nbins; iq++){
    for(int iw = 0; iw<nbins; iw++){
      myq = addhalfstep + step * (double)iq;
      myw = addhalfstep + step * (double)iw;
      int myBin_qw = susaq0q3Hist->FindBin(myq,myw);
      
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

        double w00 = (tensor->tt(myw,myq)).real();
        double w03 = (tensor->tz(myw,myq)).real();
        double w11 = (tensor->xx(myw,myq)).real();
        double w12 = -1.0*(tensor->xy(myw,myq)).imag();
        double w33 = (tensor->zz(myw,myq)).real();

        susaW00Hist->SetBinContent(myBin_qw, w00);
        susaW03Hist->SetBinContent(myBin_qw, w03);
        susaW11Hist->SetBinContent(myBin_qw, w11);
        susaW12Hist->SetBinContent(myBin_qw, w12);
        susaW33Hist->SetBinContent(myBin_qw, w33);
      }

      // Apply given overall scaling factor
      xsec *= fXSecScale;

      // Apply Jacobien (make differential)
      //xsec *= genie::utils::mec::J(myw,myq,Enu,Ml);
      //xsec *= genie::utils::mec::J(step,step,Enu,Ml);

      // output in something other than genie units.
      xsec = xsec / (1.0E-38 * units::cm2);  

      //SD 04/12/17 - apply SuSAv2 correction factor (cos^2theta_c):
      //corr_fact=0.9471182*1000.;    

      //Find bin and fill:
      susaq0q3Hist->SetBinContent(myBin_qw, xsec);

      int myBin_CthTl = susaCthTlHist->FindBin(costhl,Tl);
      susaCthTlHist->SetBinContent(myBin_CthTl, xsec);
      //susaq0q3Hist->Fill(myq,myw,xsec);

    }
  }
  susaq0q3Hist->Write();
  susaCthTlHist->Write();
  susaW00Hist->Write();
  susaW03Hist->Write();
  susaW11Hist->Write();
  susaW12Hist->Write();
  susaW33Hist->Write();
  //susaq0q3Hist->SaveAs("susaq0q3.png");
  susaq0q3->Close();    
}
//_________________________________________________________________________
void SuSAv2MECPXSec::Scanq0q3_np(void)
{
  std::cout << "Calling SuSAv2MECPXSec::Scanq0q3_np" << std::endl;


  // Get the hadron tensor pool
  HadronTensorPool& htp = HadronTensorPool::Instance();

  // Runs after loading the tensors.  But all this is hard coded, so you need a private build.
  // Can gmkspl | grep dsdq0q3 > output.dat
  
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
  HadronTensorType_t tensor_type = kHT_MEC_Fullpn;

  int nbins = 500;
  double upper = 2.0;
  double step = upper/(double)nbins;
  double addhalfstep = step/2.0;

  double corr_fact=0.0;
  double myq = 0.0;
  double myw = 0.0;

  double Q_value = 0.0;

  TFile* susaq0q3 = new TFile("susaq0q3_np.root","RECREATE");
  TH2D* susaq0q3Hist = new TH2D("SuSAv2q0q3", "SuSAv2q0q3", nbins, 0, upper, nbins, 0, upper);
  TH2D* susaCthTlHist = new TH2D("susaCthTlHist", "susaCthTlHist", 100, -1, 1, 100, 0, upper);

  TH2D*  susaW00Hist = new TH2D("susaW00Hist", "susaW00Hist", nbins, 0, upper, nbins, 0, upper); //tt
  TH2D*  susaW03Hist = new TH2D("susaW03Hist", "susaW03Hist", nbins, 0, upper, nbins, 0, upper); //tz
  TH2D*  susaW11Hist = new TH2D("susaW11Hist", "susaW11Hist", nbins, 0, upper, nbins, 0, upper); //xx
  TH2D*  susaW12Hist = new TH2D("susaW12Hist", "susaW12Hist", nbins, 0, upper, nbins, 0, upper); //xy
  TH2D*  susaW33Hist = new TH2D("susaW33Hist", "susaW33Hist", nbins, 0, upper, nbins, 0, upper); //zz

  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>( htp.GetTensor(target_pdg,
    tensor_type, "SuSAv2_MEC") );

  if(!tensor){
    std::cout << "ERROR!!!! Tensor is NULL" << std::endl;
  }

  double Q0min = tensor->q0Min();
  double Q0max = tensor->q0Max();
  double Q3min = tensor->qMagMin();
  double Q3max = tensor->qMagMax();

  std::cout << "Scanq0q3: q0 min / max is " << Q0min << " / " << Q0max << std::endl;
  std::cout << "Scanq0q3: q3 min / max is " << Q3min << " / " << Q3max << std::endl;

  std::cout << "W00 for q0=q3=100MeV  is: " << tensor->tt(0.1, 0.1) << std::endl;
  std::cout << "W00 for q0=q3=300MeV  is: " << tensor->tt(0.3, 0.3) << std::endl;
  std::cout << "W00 for q0=q3=500MeV  is: " << tensor->tt(0.5, 0.5) << std::endl;
  std::cout << "W00 for q0=q3=700MeV  is: " << tensor->tt(0.7, 0.7) << std::endl;
  std::cout << "W00 for q0=q3=900MeV  is: " << tensor->tt(0.9, 0.9) << std::endl;
  std::cout << "W00 for q0=q3=1100MeV is: " << tensor->tt(1.1, 1.1) << std::endl;
  std::cout << "W00 for q0=q3=1300MeV is: " << tensor->tt(1.3, 1.3) << std::endl << std::endl;

  std::cout << "W11 for q0=q3=100MeV  is: " << tensor->xx(0.1, 0.1) << std::endl;
  std::cout << "W11 for q0=q3=300MeV  is: " << tensor->xx(0.3, 0.3) << std::endl;
  std::cout << "W11 for q0=q3=500MeV  is: " << tensor->xx(0.5, 0.5) << std::endl;
  std::cout << "W11 for q0=q3=700MeV  is: " << tensor->xx(0.7, 0.7) << std::endl;
  std::cout << "W11 for q0=q3=900MeV  is: " << tensor->xx(0.9, 0.9) << std::endl;
  std::cout << "W11 for q0=q3=1100MeV is: " << tensor->xx(1.1, 1.1) << std::endl;
  std::cout << "W11 for q0=q3=1300MeV is: " << tensor->xx(1.3, 1.3) << std::endl << std::endl;

  std::cout << "W11 for q0=500MeV, q3=550MeV  is: " << tensor->xx(0.5, 0.55) << std::endl;
  std::cout << "W11 for q0=550MeV, q3=500MeV  is: " << tensor->xx(0.55, 0.5) << std::endl << std::endl;


  // Fixed point debugging
  double xsec_t = 0;
  double xsec_r = 0;
  double cont_t = 0;
  double w00 = 0;
  double w03 = 0;
  double w11 = 0;
  double w12 = 0;
  double w33 = 0;
  std::cout << "Calculate some fixed point xsec: " << std::endl;
  std::cout << "Enu, pmu, Tl, costhl, w, q: xsec, xsec_r, contraction" << std::endl;
  std::cout << "       HTele: w00, w03, w11, w12, w33" << std::endl;
  Enu = 1.0;
  pmu = 0.6;
  costhl = 0.8;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  //cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;


  Enu = 0.75;
  costhl = 0.5;
  pmu = 0.6;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.75;
  costhl = 0.5;
  pmu = 0.35;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 3.0;
  costhl = 0.95;
  pmu = 2.4;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 1.5;
  costhl = 0.9;
  pmu = 1.05;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.65;
  costhl = 0.902;
  pmu = 0.21065;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;


  Enu = 0.65;
  costhl = 0.902;
  pmu = 0.4;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.6;
  costhl = 0.99;
  pmu = 0.18;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.6;
  costhl = 0.99;
  pmu = 0.35;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.6;
  costhl = 0.99;
  pmu = 0.28;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t*1E+38/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r*1E+38/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 3.0; //final Enu to use


  for(int iq = 0; iq<nbins; iq++){
    for(int iw = 0; iw<nbins; iw++){
      myq = addhalfstep + step * (double)iq;
      myw = addhalfstep + step * (double)iw;
      int myBin_qw = susaq0q3Hist->FindBin(myq,myw);
      
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

        double w00 = (tensor->tt(myw,myq)).real();
        double w03 = (tensor->tz(myw,myq)).real();
        double w11 = (tensor->xx(myw,myq)).real();
        double w12 = -1.0*(tensor->xy(myw,myq)).imag();
        double w33 = (tensor->zz(myw,myq)).real();

        susaW00Hist->SetBinContent(myBin_qw, w00);
        susaW03Hist->SetBinContent(myBin_qw, w03);
        susaW11Hist->SetBinContent(myBin_qw, w11);
        susaW12Hist->SetBinContent(myBin_qw, w12);
        susaW33Hist->SetBinContent(myBin_qw, w33);
      }

      // Apply given overall scaling factor
      xsec *= fXSecScale;

      // Apply Jacobien (make differential)
      //xsec *= genie::utils::mec::J(myw,myq,Enu,Ml);
      //xsec *= genie::utils::mec::J(step,step,Enu,Ml);

      // output in something other than genie units.
      xsec = xsec / (1.0E-38 * units::cm2);  

      //SD 04/12/17 - apply SuSAv2 correction factor (cos^2theta_c):
      //corr_fact=0.9471182*1000.;    

      //Find bin and fill:
      susaq0q3Hist->SetBinContent(myBin_qw, xsec);

      int myBin_CthTl = susaCthTlHist->FindBin(costhl,Tl);
      susaCthTlHist->SetBinContent(myBin_CthTl, xsec);
      //susaq0q3Hist->Fill(myq,myw,xsec);

    }
  }
  susaq0q3Hist->Write();
  susaCthTlHist->Write();
  susaW00Hist->Write();
  susaW03Hist->Write();
  susaW11Hist->Write();
  susaW12Hist->Write();
  susaW33Hist->Write();
  //susaq0q3Hist->SaveAs("susaq0q3.png");
  susaq0q3->Close();    
}
//_________________________________________________________________________
void SuSAv2MECPXSec::Scanq0q3_electron(void)
{
  std::cout << "Calling SuSAv2MECPXSec::Scanq0q3_electron" << std::endl;


  // Get the hadron tensor pool
  HadronTensorPool& htp = HadronTensorPool::Instance();

  // Runs after loading the tensors.  But all this is hard coded, so you need a private build.
  // Can gmkspl | grep dsdq0q3 > output.dat
  
  int nu_pdg = 11;
  int target_pdg = kPdgTgtC12;
  double Enu = 3.0;
  double m_probe = PDGLibrary::Instance()->Find(11)->Mass();
  double Ml = PDGLibrary::Instance()->Find(11)->Mass();
  //double Ml = kElectronMass;
  double Tl = 0.0;
  double pmu = 0.0;
  double costhl = 0.0;
  int tensor_pdg = kPdgTgtC12;
  HadronTensorType_t tensor_type = kHT_MEC_EM;

  int nbins = 500;
  double upper = 2.0;
  double step = upper/(double)nbins;
  double addhalfstep = step/2.0;

  double corr_fact=0.0;
  double myq = 0.0;
  double myw = 0.0;

  double Q_value = -0.02;

  TFile* susaq0q3 = new TFile("susaq0q3_electron.root","RECREATE");
  TH2D* susaq0q3Hist = new TH2D("SuSAv2q0q3", "SuSAv2q0q3", nbins, 0, upper, nbins, 0, upper);
  TH2D* susaCthTlHist = new TH2D("susaCthTlHist", "susaCthTlHist", 100, -1, 1, 100, 0, upper);
  TH2D* susaCthwHist = new TH2D("susaCthwHist", "susaCthwHist", 100, -1, 1, nbins, 0, upper);

  TH2D*  susaW00Hist = new TH2D("susaW00Hist", "susaW00Hist", nbins, 0, upper, nbins, 0, upper); //tt
  TH2D*  susaW03Hist = new TH2D("susaW03Hist", "susaW03Hist", nbins, 0, upper, nbins, 0, upper); //tz
  TH2D*  susaW11Hist = new TH2D("susaW11Hist", "susaW11Hist", nbins, 0, upper, nbins, 0, upper); //xx
  TH2D*  susaW12Hist = new TH2D("susaW12Hist", "susaW12Hist", nbins, 0, upper, nbins, 0, upper); //xy
  TH2D*  susaW33Hist = new TH2D("susaW33Hist", "susaW33Hist", nbins, 0, upper, nbins, 0, upper); //zz

  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>( htp.GetTensor(target_pdg,
    tensor_type, "SuSAv2_MEC") );
  

  double Q0min = tensor->q0Min();
  double Q0max = tensor->q0Max();
  double Q3min = tensor->qMagMin();
  double Q3max = tensor->qMagMax();

  std::cout << "Scanq0q3: q0 min / max is " << Q0min << " / " << Q0max << std::endl;
  std::cout << "Scanq0q3: q3 min / max is " << Q3min << " / " << Q3max << std::endl;

  std::cout << "W00 for q0=q3=100MeV  is: " << tensor->tt(0.1, 0.1) << std::endl;
  std::cout << "W00 for q0=q3=300MeV  is: " << tensor->tt(0.3, 0.3) << std::endl;
  std::cout << "W00 for q0=q3=500MeV  is: " << tensor->tt(0.5, 0.5) << std::endl;
  std::cout << "W00 for q0=q3=700MeV  is: " << tensor->tt(0.7, 0.7) << std::endl;
  std::cout << "W00 for q0=q3=900MeV  is: " << tensor->tt(0.9, 0.9) << std::endl;
  std::cout << "W00 for q0=q3=1100MeV is: " << tensor->tt(1.1, 1.1) << std::endl;
  std::cout << "W00 for q0=q3=1300MeV is: " << tensor->tt(1.3, 1.3) << std::endl << std::endl;

  std::cout << "W11 for q0=q3=100MeV  is: " << tensor->xx(0.1, 0.1) << std::endl;
  std::cout << "W11 for q0=q3=300MeV  is: " << tensor->xx(0.3, 0.3) << std::endl;
  std::cout << "W11 for q0=q3=500MeV  is: " << tensor->xx(0.5, 0.5) << std::endl;
  std::cout << "W11 for q0=q3=700MeV  is: " << tensor->xx(0.7, 0.7) << std::endl;
  std::cout << "W11 for q0=q3=900MeV  is: " << tensor->xx(0.9, 0.9) << std::endl;
  std::cout << "W11 for q0=q3=1100MeV is: " << tensor->xx(1.1, 1.1) << std::endl;
  std::cout << "W11 for q0=q3=1300MeV is: " << tensor->xx(1.3, 1.3) << std::endl << std::endl;

  std::cout << "W11 for q0=500MeV, q3=550MeV  is: " << tensor->xx(0.5, 0.55) << std::endl;
  std::cout << "W11 for q0=550MeV, q3=500MeV  is: " << tensor->xx(0.55, 0.5) << std::endl << std::endl;


  // Fixed point debugging
  double xsec_t = 0;
  double xsec_r = 0;
  double cont_t = 0;
  double w00 = 0;
  double w03 = 0;
  double w11 = 0;
  double w12 = 0;
  double w33 = 0;
  std::cout << "Calculate some fixed point xsec: " << std::endl;
  std::cout << "Ee_probe, pe, Te, costhe, w, q: xsec, xsec_r, contraction" << std::endl;
  std::cout << "       HTele: w00, w03, w11, w12, w33" << std::endl;
  Enu = 1.0;
  pmu = 0.6;
  costhl = 0.8;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r/units::cm2;  
  //cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;


  Enu = 0.75;
  costhl = 0.5;
  pmu = 0.6;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");
  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.75;
  costhl = 0.5;
  pmu = 0.35;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 3.0;
  costhl = 0.95;
  pmu = 2.4;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 1.5;
  costhl = 0.9;
  pmu = 1.05;
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;

  Enu = 0.68;
  costhl = 0.809;
  pmu = 0.378; // w=0.3
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;


  Enu = 0.68;
  costhl = 0.809;
  pmu = 0.277; // w =0.4
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;


  Enu = 0.68;
  costhl = 0.809;
  pmu = 0.477; // w =0.2
  Tl = sqrt(pmu*pmu + Ml*Ml) - Ml;
  genie::utils::mec::Getq0q3FromTlCostl(Tl, costhl, Enu, Ml, myw, myq);
  xsec_t = tensor->dSigma_dT_dCosTheta(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_t = xsec_t/units::cm2;  
  xsec_r = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  xsec_r = xsec_r/units::cm2;  
  cont_t = tensor->contraction(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");  w00 = (tensor->tt(myw,myq)).real();
  w03 = (tensor->tz(myw,myq)).real();
  w11 = (tensor->xx(myw,myq)).real();
  w12 = -1.0*(tensor->xy(myw,myq)).imag();
  w33 = (tensor->zz(myw,myq)).real();
  std::cout << Enu << ", " << pmu << ", " << Tl << ", " << costhl << ", " << myw << ", " << myq << ": " << xsec_t << ", " << xsec_r << ", " << cont_t << std::endl;
  std::cout << "       HTele: " << w00 << ", " << w03 << ", " << w11 << ", " << w12 << ", " << w33 <<  std::endl;



  Enu = 3.0; //final Enu to use


  for(int iq = 0; iq<nbins; iq++){
    for(int iw = 0; iw<nbins; iw++){
      myq = addhalfstep + step * (double)iq;
      myw = addhalfstep + step * (double)iw;
      int myBin_qw = susaq0q3Hist->FindBin(myq,myw);
      
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

        double w00 = (tensor->tt(myw,myq)).real();
        double w03 = (tensor->tz(myw,myq)).real();
        double w11 = (tensor->xx(myw,myq)).real();
        double w12 = -1.0*(tensor->xy(myw,myq)).imag();
        double w33 = (tensor->zz(myw,myq)).real();

        susaW00Hist->SetBinContent(myBin_qw, w00);
        susaW03Hist->SetBinContent(myBin_qw, w03);
        susaW11Hist->SetBinContent(myBin_qw, w11);
        susaW12Hist->SetBinContent(myBin_qw, w12);
        susaW33Hist->SetBinContent(myBin_qw, w33);
      }

      // Apply given overall scaling factor
      xsec *= fXSecScale;

      // Apply Jacobien (make differential)
      //xsec *= genie::utils::mec::J(myw,myq,Enu,Ml);
      //xsec *= genie::utils::mec::J(step,step,Enu,Ml);

      // output in something other than genie units.
      xsec = xsec / (1.0E-38 * units::cm2);  

      //SD 04/12/17 - apply SuSAv2 correction factor (cos^2theta_c):
      //corr_fact=0.9471182*1000.;    

      //Find bin and fill:
      susaq0q3Hist->SetBinContent(myBin_qw, xsec);

      int myBin_CthTl = susaCthTlHist->FindBin(costhl,Tl);
      susaCthTlHist->SetBinContent(myBin_CthTl, xsec);
      //susaq0q3Hist->Fill(myq,myw,xsec);

      int myBin_Cthw = susaCthwHist->FindBin(costhl,myw);
      susaCthwHist->SetBinContent(myBin_Cthw, xsec);

    }
  }
  susaq0q3Hist->Write();
  susaCthTlHist->Write();
  susaCthwHist->Write();
  susaW00Hist->Write();
  susaW03Hist->Write();
  susaW11Hist->Write();
  susaW12Hist->Write();
  susaW33Hist->Write();
  susaq0q3Hist->SaveAs("susaq0q3.png");
  susaq0q3->Close();    

  std::cout << "Wait 5 seconds before continuing ..." << std::endl;
  sleep(5);

}
//_________________________________
void SuSAv2MECPXSec::buildRWhistos(void)
{
  std::cout << "Calling SuSAv2MECPXSec::buildRWhistos" << std::endl;


  // Get the hadron tensor pool
  HadronTensorPool& htp = HadronTensorPool::Instance();

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
  HadronTensorType_t tensor_type_np = kHT_MEC_Fullpn;

  int nbins = 200;
  double upper = 2.0;
  double step = upper/(double)nbins;
  double addhalfstep = step/2.0;

  double upperEnu = 5.0;
  double stepEnu = upperEnu/(double)nbins;
  double addhalfstepEnu = stepEnu/2.0;

  double corr_fact=0.0;
  double myq = 0.0;
  double myw = 0.0;
  double myenu = 0.0;

  double Q_value = 0.0;

  TFile* susaq0q3 = new TFile("SuSAMECRWHists.root","RECREATE");
  TH3D* susaEnuq0q3Hist = new TH3D("SuSAv2Elq0q3", "SuSAv2Elq0q3", nbins, 0, upperEnu, nbins, 0, upper, nbins, 0, upper);
  TH3D* susaEnuCthTlHist = new TH3D("SuSAv2ElCthTlHist", "SuSAv2ElCthTlHist", nbins, 0, upperEnu, 100, -1, 1, 100, 0, upper);

  TH3D* susaEnuq0q3Hist_pn = new TH3D("SuSAv2Elq0q3_pn", "SuSAv2Elq0q3_pn", nbins, 0, upperEnu, nbins, 0, upper, nbins, 0, upper);
  TH3D* susaEnuCthTlHist_pn = new TH3D("SuSAv2ElCthTlHist_pn", "SuSAv2ElCthTlHist_pn", nbins, 0, upperEnu, 100, -1, 1, 100, 0, upper);


  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>( htp.GetTensor(target_pdg,
    tensor_type, "SuSAv2_MEC") );

  const ValenciaHadronTensorI* tensor_pn
    = dynamic_cast<const ValenciaHadronTensorI*>( htp.GetTensor(target_pdg,
    tensor_type_np, "SuSAv2_MEC") );

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
        int myBin_qw = susaEnuq0q3Hist->FindBin(myenu,myq,myw);
        Enu = myenu;

        if(ienu%50 && iq==0 && iw==0) std::cout << "buildRWhistos: energy is " << myenu <<  ", max is " << upperEnu << std::endl;
        
        genie::utils::mec::GetTlCostlFromq0q3(myw, myq, Enu, Ml, Tl, costhl);

        double xsec = 0.0;
        double xsec_pn = 0.0;
        if(myq > Q3max || myq < Q3min || myw > Q0max || myw < Q0min){ 
          xsec = 0.0;
          xsec_pn = 0.0;
        }
        else if (Tl < 0.0){
          xsec = 0.0;
          xsec_pn = 0.0;
        }
        else { 
          xsec = tensor->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
          xsec_pn = tensor_pn->dSigma_dT_dCosTheta_rosenbluth(nu_pdg, Enu, m_probe, Tl, costhl, Ml, Q_value);
        }

        // Apply given overall scaling factor
        xsec *= fXSecScale;
        xsec_pn *= fXSecScale;

        // output in something other than genie units.
        xsec = xsec / (1.0E-38 * units::cm2);  
        xsec_pn = xsec_pn / (1.0E-38 * units::cm2);  

        //Find bin and fill:
        susaEnuq0q3Hist->SetBinContent(myBin_qw, xsec);
        susaEnuq0q3Hist_pn->SetBinContent(myBin_qw, xsec_pn);

        int myBin_CthTl = susaEnuCthTlHist->FindBin(myenu, costhl,Tl);
        susaEnuCthTlHist->SetBinContent(myBin_CthTl, xsec);
        susaEnuCthTlHist_pn->SetBinContent(myBin_CthTl, xsec_pn);
        //susaq0q3Hist->Fill(myq,myw,xsec);

      }
    }
  }

  susaEnuq0q3Hist->Write();
  susaEnuCthTlHist->Write();

  susaEnuq0q3Hist_pn->Write();
  susaEnuCthTlHist_pn->Write();

  //susaq0q3Hist->SaveAs("susaq0q3.png");
  susaq0q3->Close();    
}

