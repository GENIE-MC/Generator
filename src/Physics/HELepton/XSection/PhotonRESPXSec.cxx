//____________________________________________________________________________
/*
*/
//____________________________________________________________________________

#include "Physics/HELepton/XSection/PhotonStrucFunc.h"
#include "Physics/HELepton/XSection/PhotonRESPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
PhotonRESPXSec::PhotonRESPXSec() :
XSecAlgorithmI("genie::PhotonRESPXSec")
{
  born = new Born();
}
//____________________________________________________________________________
PhotonRESPXSec::PhotonRESPXSec(string config) :
XSecAlgorithmI("genie::PhotonRESPXSec", config)
{

}
//____________________________________________________________________________
PhotonRESPXSec::~PhotonRESPXSec()
{

}
//____________________________________________________________________________
double PhotonRESPXSec::XSec(
   const Interaction * interaction, KinePhaseSpace_t kps) const
{

  if(! this -> ValidProcess    (interaction) ) return 0.;

  // Load SF tables 
  PhotonStrucFunc * sf_tbl = PhotonStrucFunc::Instance();

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const XclsTag &      xclstag    = interaction -> ExclTag();

  int probepdg = init_state.ProbePdg();
  int loutpdg  = xclstag.FinalLeptonPdg();
  int tgtpdg   = init_state.Tgt().HitNucPdg();

  double mlin  = 0.;
  if      (pdg::IsNuE  (TMath::Abs(probepdg))) mlin = kElectronMass;
  else if (pdg::IsNuMu (TMath::Abs(probepdg))) mlin = kMuonMass;
  else if (pdg::IsNuTau(TMath::Abs(probepdg))) mlin = kTauMass;
  double mlout = interaction->FSPrimLepton()->Mass(); //mass of charged lepton

  double Mnuc = init_state.Tgt().HitNucMass();

  double E = init_state.ProbeE(kRfLab);
  double s = (2 * Mnuc * E + Mnuc*Mnuc);

  double n1 = kinematics.GetKV(kKVn1);
  double n2 = kinematics.GetKV(kKVn2);

  double xmin = fQ2PDFmin/2./E/Mnuc;
  double x = TMath::Exp( TMath::Log(xmin) + (TMath::Log(1.0)-TMath::Log(xmin))*n2 );

  if (x<fxPDFmin) return 0.;

  double s_r = x*s;
  double t_r = born->GetT(0.,mlin,interaction->FSPrimLepton()->Mass(),0.,s_r,n1);

  double xsec = kPi/4./(s_r-Mnuc*Mnuc) * sf_tbl->EvalSF(tgtpdg,probepdg,x) * (TMath::Log(1.0)-TMath::Log(xmin)) ;
  
  if ( pdg::IsPion(loutpdg) ) {
    if ( TMath::Sqrt(s_r)<fWmin ) return 0.; // The W limit is because hadronization might have issues at low W (as in PYTHIA6).
    xsec *= 64.41/10.63;    
  }

  double ME = 0.;
  if ( TMath::Abs(loutpdg)+1 == TMath::Abs(probepdg) ) ME = PXSecCCRNC(s,t,mlin*mlin,mlout*mlout);
  else                                                 ME = PXSecCCR  (s,t,mlin*mlin,mlout*mlout); 
  xsec *= TMath::Max(0.,ME);
   
  if(kps!=kPSn1n2fE) {
      LOG("PhotonRESPXSec", pWARN)
          << "Doesn't support transformation from "
          << KinePhaseSpace::AsString(kPSn1n2fE) << " to "
          << KinePhaseSpace::AsString(kps);
      xsec = 0;
  }

  // If requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Compute nuclear cross section (simple scaling here, corrections must have been included in the structure functions)
  int NNucl = (pdg::IsProton(tgtpdg)) ? init_state.Tgt().Z() : init_state.Tgt().N(); 
  xsec *= NNucl; 

  LOG("PhotonRESPXSec", pINFO) << "dxsec/dn1dn2 (E= " << E << ", n1= " << n1 << ", n2=" << n2 << ") = " << xsec;

  return xsec;

}
//____________________________________________________________________________
double PhotonRESPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool PhotonRESPXSec::ValidProcess(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsPhotonRES()) return false;

  const InitialState & init_state = interaction -> InitState();
  if(!pdg::IsLepton(init_state.ProbePdg())) return false;

  if(! init_state.Tgt().HitNucIsSet()) return false;

  int hitnuc_pdg = init_state.Tgt().HitNucPdg();
  if(!pdg::IsNeutronOrProton(hitnuc_pdg)) return false;
  
  return true;
}
//____________________________________________________________________________


//____________________________________________________________________________
void PhotonRESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonRESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonRESPXSec::LoadConfig(void)
{

  //-- load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  GetParam( "Xsec-Wmin", fWmin ) ;

  GetParam("Q2Grid-Min", fQ2PDFmin );
  GetParam("XGrid-Min",  fxPDFmin );

}