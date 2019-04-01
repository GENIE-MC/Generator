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
#include "Physics/QuasiElastic/XSection/SuSAv2QEXSec.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
// For q3q0 scan verification and RW histos:
#include <TH1F.h>
#include <TH3D.h>
#include <TFile.h>
#include "Framework/Conventions/Constants.h"

using namespace genie;

//_________________________________________________________________________
SuSAv2QEXSec::SuSAv2QEXSec() : XSecAlgorithmI("genie::SuSAv2QEXSec")
{
}
//_________________________________________________________________________
SuSAv2QEXSec::SuSAv2QEXSec(string config)
  : XSecAlgorithmI("genie::SuSAv2QEXSec", config)
{
}
//_________________________________________________________________________
SuSAv2QEXSec::~SuSAv2QEXSec()
{
}
//_________________________________________________________________________
double SuSAv2QEXSec::XSec(const Interaction* interaction,
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
    tensor_type = kHT_QE_Full;
  }
  else {
    // If the probe is not a neutrino, assume that it's an electron
    tensor_type = kHT_QE_EM;
  }

  // The SuSAv2-1p1h hadron tensors are defined using the same conventions
  // as the Valencia MEC (and SuSAv2-MEC) model, so we can use the same sort of tensor
  // object to describe them.
  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>( htp.GetTensor(target_pdg,
    tensor_type, fHadronTensorTableName) );

  /// \todo Add scaling of the hadron tensor elements if an exact match
  /// for the requested target PDG code cannot be found in the hadron tensor
  /// table for this model.
  /// Also add the different pair configurations for e-scattering

  // If retrieving the tensor failed, complain and return zero
  if ( !tensor ) {
    LOG("SuSAv2QE", pWARN) << "Failed to load a hadronic tensor for the"
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
  // SD: The Q-Value essentially corrects q0 to account for nuclear
  // binding energy in the Valencia model but this effect is already
  // in Guille's tensors so I'll set it to 0.

  int nu_pdg = interaction->InitState().ProbePdg();
  double Q_value = 0.;


  // Compute the cross section using the hadron tensor
  double xsec = tensor->dSigma_dT_dCosTheta_rosenbluth(interaction, Q_value);
  LOG("SuSAv2QE", pDEBUG) << "XSec in cm2 / neutron is  " << xsec/(units::cm2);

  // Apply the squared CKM matrix element Vud as a correction factor
  // (not included in the tabulated tensor element values)
  // NOT NEEDED IF USING THE ROSENBLUTH FORMALISM
  //xsec *= fVud2;

  // Currently the hadron tensors are per neutron, but the calculation above
  // assumes they are per atom. Need to adjust for this
  xsec *= interaction->InitState().Tgt().Z();
  LOG("SuSAv2QE", pDEBUG) << "XSec in cm2 / atom is  " << xsec/(units::cm2);

  // Apply given overall scaling factor
  xsec *= fXSecScale;

  if ( kps != kPSTlctl ) {
    LOG("SuSAv2QE", pWARN)
      << "Doesn't support transformation from "
      << KinePhaseSpace::AsString(kPSTlctl) << " to "
      << KinePhaseSpace::AsString(kps);
    xsec = 0.;
  }

  return xsec;
}
//_________________________________________________________________________
double SuSAv2QEXSec::Integral(const Interaction* interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//_________________________________________________________________________
bool SuSAv2QEXSec::ValidProcess(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if(!proc_info.IsQuasiElastic()) return false;

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  bool isP   = pdg::IsProton(nuc);
  bool isN   = pdg::IsNeutron(nuc);
  bool isnu  = pdg::IsNeutrino(nu);
  bool isnub = pdg::IsAntiNeutrino(nu);

  bool prcok = proc_info.IsWeakCC() && ((isP&&isnub) || (isN&&isnu));
  if(!prcok) return false;

  return true;
}
//_________________________________________________________________________
void SuSAv2QEXSec::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SuSAv2QEXSec::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//_________________________________________________________________________
void SuSAv2QEXSec::LoadConfig(void)
{
  // CKM matrix element connecting the up and down quarks
  // (not included in the tabulated SuSAv2-QE hadron tensors)
  double Vud;
  GetParam("CKM-Vud", Vud);
  fVud2 = std::pow(Vud, 2);

  // Cross section scaling factor
  GetParam( "QEL-CC-XSecScale", fXSecScale ) ;

  // Name of the hadron tensor table to use for this model
  GetParamDef("HadronTensorTableName", fHadronTensorTableName,
    std::string("SuSAv2_QE"));

   // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  //uncomment this to make histo output on a 2D grid for validation
  this->Scanq0q3();
  //this->Scanq0q3_np();
  //this->Scanq0q3_electron();

  //uncomment this to make 3D histograms to be used for RW-ing
  //this->buildRWhistos();
}
//_________________________________________________________________________
void SuSAv2QEXSec::Scanq0q3(void)
{
  std::cout << "Calling SuSAv2QEXSec::Scanq0q3" << std::endl;


  // Get the hadron tensor pool
  HadronTensorPool& htp = HadronTensorPool::Instance();

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
  HadronTensorType_t tensor_type = kHT_QE_Full;

  int nbins = 500;
  double upper = 2.0;
  double step = upper/(double)nbins;
  double addhalfstep = step/2.0;

  double corr_fact=0.0;
  double myq = 0.0;
  double myw = 0.0;

  double Q_value = 0.0;

  TFile* susaq0q3 = new TFile("susaq0q3_QE.root","RECREATE");
  TH2D* susaq0q3Hist = new TH2D("SuSAv2q0q3", "SuSAv2q0q3", nbins, 0, upper, nbins, 0, upper);
  TH2D* susaCthTlHist = new TH2D("susaCthTlHist", "susaCthTlHist", 100, -1, 1, 100, 0, upper);

  TH2D*  susaW00Hist = new TH2D("susaW00Hist", "susaW00Hist", nbins, 0, upper, nbins, 0, upper); //tt
  TH2D*  susaW03Hist = new TH2D("susaW03Hist", "susaW03Hist", nbins, 0, upper, nbins, 0, upper); //tz
  TH2D*  susaW11Hist = new TH2D("susaW11Hist", "susaW11Hist", nbins, 0, upper, nbins, 0, upper); //xx
  TH2D*  susaW12Hist = new TH2D("susaW12Hist", "susaW12Hist", nbins, 0, upper, nbins, 0, upper); //xy
  TH2D*  susaW33Hist = new TH2D("susaW33Hist", "susaW33Hist", nbins, 0, upper, nbins, 0, upper); //zz

  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>( htp.GetTensor(target_pdg,
    tensor_type, "SuSAv2_QE") );

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
  //cont_t = genie::utils::mec::OldTensorContraction(nu_pdg, target_pdg, Enu, Ml, Tl, costhl, target_pdg, tensor_type, "SuSAv2_MEC");
  cont_t = 0.0; 
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
void SuSAv2QEXSec::Scanq0q3_np(void)
{
  std::cout << "Calling SuSAv2QEXSec::Scanq0q3_np" << std::endl;


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
void SuSAv2QEXSec::Scanq0q3_electron(void)
{
  std::cout << "Calling SuSAv2QEXSec::Scanq0q3_electron" << std::endl;


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
void SuSAv2QEXSec::buildRWhistos(void)
{
  std::cout << "Calling SuSAv2QEXSec::buildRWhistos" << std::endl;


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

  TFile* susaq0q3 = new TFile("SuSAMECRWHists.root","RECREATE");
  TH3D* susaEnuq0q3Hist = new TH3D("SuSAv2Elq0q3", "SuSAv2Elq0q3", nbins, 0, upperEnu, nbins, 0, upper, nbins, 0, upper);
  TH3D* susaEnuCthTlHist = new TH3D("SuSAv2ElCthTlHist", "SuSAv2ElCthTlHist", nbins, 0, upperEnu, 100, -1, 1, 100, 0, upper);



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

  std::cout << "buildRWhistos: q0 min / max is " << Q0min << " / " << Q0max << std::endl;
  std::cout << "buildRWhistos: q3 min / max is " << Q3min << " / " << Q3max << std::endl;


  for(int ienu = 0; ienu<nbins; ienu++){
    for(int iq = 0; iq<nbins; iq++){
      for(int iw = 0; iw<nbins; iw++){
        myenu = addhalfstepEnu + stepEnu * (double)ienu;
        myq = addhalfstep + step * (double)iq;
        myw = addhalfstep + step * (double)iw;
        int myBin_qw = susaEnuq0q3Hist->FindBin(myenu,myq,myw);

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
        susaEnuq0q3Hist->SetBinContent(myBin_qw, xsec);

        int myBin_CthTl = susaEnuCthTlHist->FindBin(myenu, costhl,Tl);
        susaEnuCthTlHist->SetBinContent(myBin_CthTl, xsec);
        //susaq0q3Hist->Fill(myq,myw,xsec);

      }
    }
  }

  susaEnuq0q3Hist->Write();
  susaEnuCthTlHist->Write();

  //susaq0q3Hist->SaveAs("susaq0q3.png");
  susaq0q3->Close();    
}