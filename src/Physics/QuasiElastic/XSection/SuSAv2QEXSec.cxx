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
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"

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
  int tensor_pdg = target_pdg;
  int A_request = pdg::IonPdgCodeToA(target_pdg);
  int Z_request = pdg::IonPdgCodeToZ(target_pdg);
  bool need_to_scale = false;

  HadronTensorType_t tensor_type = kHT_Undefined;
  if ( pdg::IsNeutrino(probe_pdg) || pdg::IsAntiNeutrino(probe_pdg) ) {
    tensor_type = kHT_QE_Full;
  }
  else {
    // If the probe is not a neutrino, assume that it's an electron
    tensor_type = kHT_QE_EM;
  }

  /// \todo Add more hadron tensors so this scaling is not so terrible
  // At the moment all we have is Carbon so this is all just a place holder ... 
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
      tensor_pdg = kPdgTgtC12;    
  }
  else if(A_request >= 22) {
      tensor_pdg = kPdgTgtC12;
      LOG("SuSAv2MEC", pWARN) << "Dangerous scaling from " << target_pdg << " to " << tensor_pdg;    
      LOG("SuSAv2MEC", pWARN) << "Hadron tensors for such heavy elements not ready yet";    
      LOG("SuSAv2MEC", pWARN) << "Proceed at your own peril!";    
  }

  if(tensor_pdg != target_pdg) need_to_scale = true;

  // The SuSAv2-1p1h hadron tensors are defined using the same conventions
  // as the Valencia MEC (and SuSAv2-MEC) model, so we can use the same sort of tensor
  // object to describe them.
  const ValenciaHadronTensorI* tensor
    = dynamic_cast<const ValenciaHadronTensorI*>( htp.GetTensor(tensor_pdg,
    tensor_type, fHadronTensorTableName) );


  // If retrieving the tensor failed, complain and return zero
  if ( !tensor ) {
    LOG("SuSAv2QE", pWARN) << "Failed to load a hadronic tensor for the"
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
    double scaleFact = (KF_tgt/KF_ten); // A-scaling already applied in section above
    xsec *= scaleFact;
  }


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
  // (not included in the tabulated SuSAv2-MEC hadron tensors)
  // No longer needed as using Rosenbluth formalism
  //double Vud;
  //GetParam("CKM-Vud", Vud);
  //fVud2 = std::pow(Vud, 2);

  // Cross section scaling factor
  GetParam( "QEL-CC-XSecScale", fXSecScale ) ;

  // Name of the hadron tensor table to use for this model
  GetParamDef("HadronTensorTableName", fHadronTensorTableName,
    std::string("SuSAv2_QE"));

   // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
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


//Details on what form factors are used for SuSAv2:

/*
NOTE: MAQE = 1.032 using dipole form factor


       subroutine GKex(Q2,GEPROT,GMPROT,GENEUT,GMNEUT)
       implicit real*8 (a-h,o-z)
C
C---------------------------------------------------------------
C                           All masses, etc. in GeV
       q2t=-q2 !Esto seria Q^2
              xmnucleon=0.939D0
C                    nucleon mass
              xmrho=0.776D0
C                    rho mass
C
C                           irho=0 : no rho width
C                           irho=1 : with rho width
       xirho=1.      
              xmrho1=xmrho-xirho*0.03465D0
C                    mass rho1
              xmrho2=xmrho-xirho*0.04374D0
C                    mass rho2
              alpha1=xirho*0.0781808D0
C                    parameter alpha1
              alpha2=xirho*0.0632907D0
C                    parameter alpha2
              qq1=0.3176D0
C                    parameter q1**2
              qq2=0.1422D0
C                    parameter q2**2
              xmrhop=1.45D0
C                    rho-prime mass
              xmomega=0.784D0
C                    omega mass
              xmomegap=1.419D0
C                    omega-prime mass
              xmphi=1.019D0
C                    phi mass
C
C----------------------------------------------------------------
C
              xkaps=-0.12D0
C                    isoscalar anomalous mag. mom.
              xkapv=3.706D0        
C                    isovector anomalous mag. mom.
              xkaprho=5.51564D0
C                    rho anomalous moment
              xkaprhop=12.0D0
C                    rho-prime anomalous moment
              xkapomega=0.4027D0
C                    omega anomalous moment
              xkapomegap=-2.973D0
C                    omega-prime anomalous moment
              xkapphi=0.01D0
C                    phi anomalous moment
C
C----------------------------------------------------------------
C
              fgrho=0.5596D0
C                    rho coupling
              fgrhop=0.007208D0
C                    rho-prime coupling
              fgomega=0.7021D0
C                    omega coupling
              fgomegap=0.164D0
C                    omega-prime coupling
              fgphi=-0.1711D0
C                    phi coupling
C
C----------------------------------------------------------------
C
              xlam1=0.93088D0
C                    parameter lambda1
              xlam2=2.6115D0
C                    parameter lambda2
              xlamd=1.181D0
C                    parameter lambdaD
              xlamqcd=0.15D0
C                    parameter lambdaQCD
              arat= DLOG(xlamd**2./xlamqcd**2.)
              xmuphi=0.2D0
C                    parameter mu-phi
C
C------------------------------------------------------------------
c
              tau=q2t/(4.0D0*xmnucleon**2.)
C                    q2 is 4-momentum q**2; tau is the usual
              femrho1=xmrho1**2./(xmrho1**2.+q2t)
              femrho2=xmrho2**2./(xmrho2**2.+q2t)
              femrhop=xmrhop**2./(xmrhop**2.+q2t)
              femomega=xmomega**2./(xmomega**2.+q2t)
              femomegap=xmomegap**2./(xmomegap**2.+q2t)
              femphi=xmphi**2./(xmphi**2.+q2t)
C
C      In the following the notation should be obvious: 
C      e.g., f1somega means F1, isoscalar (s), omega piece, etc.
C
c             
              q2tilda=q2t*DLOG((xlamd**2.+q2t)/xlamqcd**2.)/arat
              f1=xlam1**2./(xlam1**2.+q2tilda)
              f2=xlam2**2./(xlam2**2.+q2tilda)
              fhad1=f1*f2
              fhad2=f1*f2**2.
              fhad1s=fhad1*(q2t/(xlam1**2.+q2t))**1.5
              fmuphi=(xmuphi**2.+q2t)/xmuphi**2.
              fhad2s=fhad2*(fmuphi*xlam1**2./(xlam1**2.+q2t))**1.5
              fqcd=xlamd**2./(xlamd**2.+q2tilda)
              fhad1qcd=fqcd*f2
              fhad2qcd=fqcd*f2**2.
c
C
              zzf1s=1.-fgomega-fgomegap
       zzf2s=xkaps-xkapomega*fgomega-xkapomegap*fgomegap-xkapphi*fgphi
              zzf1v=1.-fgrho-fgrhop
              zzf2v=xkapv-xkaprho*fgrho-xkaprhop*fgrhop
c
c
c
              f1somega=fgomega*femomega*fhad1
              f1somegap=fgomegap*femomegap*fhad1
              f1sphi= fgphi*femphi*fhad1s
              f1spqcd= zzf1s*fhad1qcd
              f2somega=xkapomega*fgomega*femomega*fhad2
              f2somegap= xkapomegap*fgomegap*femomegap*fhad2
              f2sphi= xkapphi*fgphi*femphi*fhad2s
              f2spqcd= zzf2s*fhad2qcd
              width1=1.-alpha1+alpha1/(1.+q2t/qq1)**2.
              width2=1.-alpha2+alpha2/(1.+q2t/qq2)
              f1vrho=fgrho*femrho1*fhad1*width1
              f1vrhop= fgrhop*femrhop*fhad1
              f1vpqcd= zzf1v*fhad1qcd
              f2vrho=xkaprho*fgrho*femrho2*fhad2*width2
              f2vrhop= xkaprhop*fgrhop*femrhop*fhad2
              f2vpqcd= zzf2v*fhad2qcd
c
C
c      
              f1s=f1somega+f1somegap+f1sphi+f1spqcd
              f2s=f2somega+f2somegap+f2sphi+f2spqcd
              f1v=f1vrho+f1vrhop+f1vpqcd
              f2v=f2vrho+f2vrhop+f2vpqcd
c
C
c             
              f1prho=0.5D0*(f1vrho)
              f1prhop= 0.5D0*(f1vrhop)
              f1pomega= 0.5D0*(f1somega)
              f1pomegap= 0.5D0*(f1somegap)
              f1pphi= 0.5D0*(f1sphi)
              f1ppqcd= 0.5D0*(f1spqcd+f1vpqcd)
              f1p=f1prho+f1prhop+f1pomega+f1pomegap+f1pphi+f1ppqcd
C
              f1nrho=0.5D0*(-f1vrho)
              f1nrhop= 0.5D0*(-f1vrhop)
              f1nomega= 0.5D0*(f1somega)
              f1nomegap= 0.5D0*(f1somegap)
              f1nphi= 0.5D0*(f1sphi)
              f1npqcd= 0.5D0*(f1spqcd-f1vpqcd)
              f1n=f1nrho+f1nrhop+f1nomega+f1nomegap+f1nphi+f1npqcd
C
              f2prho=0.5D0*(f2vrho)
              f2prhop= 0.5D0*(f2vrhop)
              f2pomega= 0.5D0*(f2somega)
              f2pomegap= 0.5D0*(f2somegap)
              f2pphi= 0.5D0*(f2sphi)
              f2ppqcd= 0.5D0*(f2spqcd+f2vpqcd)
              f2p=f2prho+f2prhop+f2pomega+f2pomegap+f2pphi+f2ppqcd
C
              f2nrho=0.5D0*(-f2vrho)
              f2nrhop= 0.5D0*(-f2vrhop)
              f2nomega= 0.5D0*(f2somega)
              f2nomegap= 0.5D0*(f2somegap)
              f2nphi= 0.5D0*(f2sphi)
              f2npqcd= 0.5D0*(f2spqcd-f2vpqcd)
              f2n=f2nrho+f2nrhop+f2nomega+f2nomegap+f2nphi+f2npqcd
C
c             
              gevrho=f1vrho-tau*f2vrho
              gevrhop=f1vrhop-tau*f2vrhop
              gesomega=f1somega-tau*f2somega
              gesomegap=f1somegap-tau*f2somegap
              gesphi=f1sphi-tau*f2sphi
              gespqcd=f1spqcd-tau*f2spqcd
              gevpqcd=f1vpqcd-tau*f2vpqcd
              ges=gesomega+gesomegap+gesphi+gespqcd
              gev=gevrho+gevrhop+gevpqcd
C
              gmvrho=f1vrho+f2vrho
              gmvrhop=f1vrhop+f2vrhop
              gmsomega=f1somega+f2somega
              gmsomegap=f1somegap+f2somegap
              gmsphi=f1sphi+f2sphi
              gmspqcd=f1spqcd+f2spqcd
              gmvpqcd=f1vpqcd+f2vpqcd
              gms=gmsomega+gmsomegap+gmsphi+gmspqcd
              gmv=gmvrho+gmvrhop+gmvpqcd
C
              geprho= 0.5D0*(gevrho)
              geprhop= 0.5D0*(gevrhop)
              gepomega= 0.5D0*(gesomega)
              gepomegap= 0.5D0*(gesomegap)
              gepphi= 0.5D0*(gesphi)
              geppqcd= 0.5D0*(gespqcd+gevpqcd)
              xgep=geprho+geprhop+gepomega+gepomegap+gepphi+geppqcd
C
              genrho= 0.5D0*(-gevrho)
              genrhop= 0.5D0*(-gevrhop)
              genomega= 0.5D0*(gesomega)
              genomegap= 0.5D0*(gesomegap)
              genphi= 0.5D0*(gesphi)
              genpqcd= 0.5D0*(gespqcd-gevpqcd)
              xgen=genrho+genrhop+genomega+genomegap+genphi+genpqcd
C
              gmprho= 0.5D0*(gmvrho)
              gmprhop= 0.5D0*(gmvrhop)
              gmpomega= 0.5D0*(gmsomega)
              gmpomegap= 0.5D0*(gmsomegap)
              gmpphi= 0.5D0*(gmsphi)
              gmppqcd= 0.5D0*(gmspqcd+gmvpqcd)
              xgmp=gmprho+gmprhop+gmpomega+gmpomegap+gmpphi+gmppqcd
C
              gmnrho= 0.5D0*(-gmvrho)
              gmnrhop= 0.5D0*(-gmvrhop)
              gmnomega= 0.5D0*(gmsomega)
              gmnomegap= 0.5D0*(gmsomegap)
              gmnphi= 0.5D0*(gmsphi)
              gmnpqcd= 0.5D0*(gmspqcd-gmvpqcd)
              xgmn=gmnrho+gmnrhop+gmnomega+gmnomegap+gmnphi+gmnpqcd
c
C Al final se obtienen los factores de forma electricos y magneticos de Sachs para protones y neutrones, y para el caso isovector (CC neutrino)
       GEPROT=xgep
       GENEUT=xgen
       GMPROT=xgmp
       GMNEUT=xgmn
  GE1=GEPROT-GENEUT
  GM1=GMPROT-GMNEUT
              
              return
              end

*/