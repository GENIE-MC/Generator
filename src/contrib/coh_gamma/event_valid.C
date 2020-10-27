
#include "TString.h" 
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

#include "Framework/Ntuple/NtpMCTreeHeader.h" 
#include "Framework/Ntuple/NtpMCEventRecord.h" 
#include "Framework/EventGen/EventRecord.h" 
#include "Framework/ParticleData/BaryonResUtils.h"

#include "Framework/GHEP/GHepParticle.h"       
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/GHEP/GHepStatus.h"

#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"


using namespace genie ;


double LepPhi( const TLorentzVector & probe, TLorentzVector lep ) ;
double GammaPhi( const TLorentzVector & probe, TLorentzVector gamma, TLorentzVector lep ) ;
bool SelEvent( double theta_l, double theta_g, double phi_g ) ;
void GraphIntegratedCOHGammaXsec( double nevts, TH1 * h_E_g, TFile & file ) ;

/////////////
// Global constants

double g_theta_l_max = 0.6 ;
double g_theta_l_min = 0.4 ;
double g_theta_g_max = 0.6 ;
double g_theta_g_min = 0.4 ;
double g_phi_g_max   = 3.0 ;
double g_phi_g_min   = 2.8 ;
double g_Ev      = 0.5 ;
int g_target     = 1000060120 ;
// COH Gamma xsec config set ("Default" or "Consistent")
std::string config = "Default" ;


void event_valid( TString in_file_name  = "../event_gen/gntp.0.ghep.root" , 
		 TString out_file_name = "" ) {

  
  if ( out_file_name == "" ) {
    out_file_name = in_file_name ; 
    out_file_name.ReplaceAll( ".ghep.root", 
			      ".evtplots.root" ) ;
  }
    
  TFile in_file( in_file_name ) ; 
    
  NtpMCTreeHeader * header = dynamic_cast<NtpMCTreeHeader*> (in_file.Get("header"));

  // Get the GENIE GHEP tree and set its branch address
  TTree * in_tree = dynamic_cast<TTree*> (in_file.Get("gtree"));
  NtpMCEventRecord * mcrec = 0;
  in_tree->SetBranchAddress("gmcrec", & mcrec);

  
  std::map<std::string, TH1*> hists ;

  // plots we need
  TH1* h_theta_gamma = hists["theta_gamma"] = new TH1D("h_theta_gamma", "#theta_{#gamma};#theta_{#gamma} [rad]",
							80, 0., TMath::Pi() ) ; 
  TH1* h_theta_lep = hists["theta_lep"] = new TH1D("h_theta_lep", "#theta_{l};#theta_{l} [rad]",
						   80, 0., TMath::Pi() ) ; 
  TH1* h_theta_recoil = hists["theta_nuc"] = new TH1D("h_theta_recoil", "#theta_{recoil};#theta_{recoil} [rad]",
						   80, 0., TMath::Pi() ) ; 

  TH1* h_phi_gamma = hists["phi_gamma"] = new TH1D("h_phi_gamma", "#phi_{#gamma};#phi_{#gamma} [rad]",
							100, -TMath::Pi(), TMath::Pi() ) ; 
  TH1* h_phi_lep = hists["phi_lep"] = new TH1D("h_phi_lep", "#phi_{l};#phi_{#phi} [rad]",
							100, -TMath::Pi(), TMath::Pi() ) ; 

  TH1* h_E_l = hists["E_l"] = new TH1D("h_E_l", "E_{l};E_{l} [GeV]",
				       80, 0., 2. ) ; 

  TH1* h_E_g = hists["E_g"] = new TH1D("h_E_g", "E_{#gamma};E_{#gamma} [GeV]",
				       80, 0., 2. ) ; 
  TH1* h_E_nu = hists["E_nu"] = new TH1D("h_E_nu", "E_{#nu};E_{#nu} [GeV]",
        			       100, 0., 10. ) ; 


  double nevts = in_tree->GetEntries() ;
  size_t evts = 0;
  size_t sel_evts = 0;

  // Event loop
  for(Long64_t i=0; i < in_tree->GetEntries(); i++) {
 
    in_tree->GetEntry(i);
    
    EventRecord & event = *(mcrec->event);
    
    const Interaction & inter = *( event.Summary() ) ;
    
    const ProcessInfo & proc_info = inter.ProcInfo() ;
    
    if ( proc_info.IsCoherentProduction() ) {
    
      const XclsTag & xcl_tag = inter.ExclTag() ;

      if ( xcl_tag.NSingleGammas() == 1 ) {
	
	TLorentzVector probe  ( * event.Probe() -> P4() ) ;
	TLorentzVector lep    ( * event.FinalStatePrimaryLepton() -> P4() ) ;
	TLorentzVector gamma  ( * event.Particle(3) -> P4() ) ;
	TLorentzVector recoil ( * event.Particle(4) -> P4() ) ;

        double theta_l = lep.Angle( probe.Vect() ) ;
        double theta_g = gamma.Angle( probe.Vect() ) ;
        double phi_g = GammaPhi( probe, gamma, lep ) ;

	h_theta_gamma -> Fill( theta_g ) ;
	h_theta_lep   -> Fill( theta_l ) ;
	h_theta_recoil   -> Fill( recoil.Angle( probe.Vect() ) ) ;

	h_phi_gamma -> Fill( phi_g ) ;			     
	h_phi_lep -> Fill( LepPhi( probe, lep ) ) ;			     

	h_E_l -> Fill( lep.E() ) ;
	h_E_nu -> Fill( probe.E() ) ;
      
        if ( !SelEvent( theta_l, theta_g, phi_g ) ) {  mcrec->Clear(); continue; }

	h_E_g -> Fill( gamma.E() ) ;

        sel_evts += 1;
	
      } // single gamma

    } // coherent production

      mcrec->Clear() ;
  } // event loop

   
  TFile out_file ( out_file_name, "RECREATE" ) ;
  out_file.cd() ;
  
  GraphIntegratedCOHGammaXsec( nevts, h_E_g, out_file ) ;

  for ( auto & h : hists ) {
    h.second -> Write() ;
    delete h.second ;
  }

  out_file.Close() ;

  std::cout << "Processed " << nevts << " and selected " << sel_evts << " events." << std::endl;

}



double LepPhi( const TLorentzVector & probe, TLorentzVector lep ) {


  TVector3 probe_dir = probe.Vect().Unit() ;
  lep.RotateUz( probe_dir ) ;
  return lep.Phi() ;

}

double GammaPhi( const TLorentzVector & probe, TLorentzVector gamma, TLorentzVector lep ) {

  TVector3 probe_dir = probe.Vect().Unit() ;
  gamma.RotateUz( probe_dir ) ;
  lep.RotateUz( probe_dir ) ;

  return gamma.DeltaPhi( lep ) ;

}

bool SelEvent( double theta_l, double theta_g, double phi_g ) {

  if ( ( theta_l > g_theta_l_min && theta_l < g_theta_l_max ) &&
       ( theta_g > g_theta_g_min && theta_g < g_theta_g_max ) &&
       ( phi_g   > g_phi_g_min && phi_g     < g_phi_g_max ) ) return true;

  return false;
  
}

void GraphIntegratedCOHGammaXsec( double nevts, TH1 * h_E_g, TFile & file ) {
  
  // -- request the COH NC Gamma mode 
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("AlvarezRusoSalaCOHGammaPXSec",config) ;
  const Algorithm * algXsec = algf->GetAlgorithm(id);
  const XSecAlgorithmI * fXSec = dynamic_cast<const XSecAlgorithmI *>(algXsec);
   
  Interaction * in = Interaction::COHNC(g_target, 14, 22, g_Ev);
  utils::gsl::d4Xsec_dEgdThetaldThetagdPhig func = utils::gsl::d4Xsec_dEgdThetaldThetagdPhig( fXSec, in ) ;

  double theta_l = ( g_theta_l_max + g_theta_l_min ) / 2. ;
  double theta_g = ( g_theta_g_max + g_theta_g_min ) / 2. ;
  double phi_g   = ( g_phi_g_max + g_phi_g_min ) / 2. ;

  double wtheta_l = ( g_theta_l_max - g_theta_l_min ) ;
  double wtheta_g = ( g_theta_g_max - g_theta_g_min ) ;
  double wphi_g   = ( g_phi_g_max - g_phi_g_min ) ;

  int nbin = h_E_g -> GetNbinsX() ;
  double xsec_i_arr[nbin];
  double xsec_arr[nbin];
  double eg_arr[nbin];

  for ( int b = 1; b < nbin; b++ ) {
    double Eg = h_E_g -> GetXaxis() -> GetBinCenter( b ) ;
    double Ni = h_E_g -> GetBinContent( b ) ;
    double wbin = h_E_g -> GetBinWidth( b ) ;

    eg_arr[b-1] = Eg ;
    std::array<double, 4> point = { Eg, theta_l, theta_g, phi_g } ; 
    xsec_i_arr[b-1] =  1.2301927672e-14 * ( Ni / nevts ) / ( wbin * wtheta_l * wtheta_g * wphi_g ) ; // sigma(Ev) from spine file for 1GeV
    xsec_arr[b-1] = func( point.data() ) ;
    std::cout << "Eg = " << Eg << "Sigma functor = " << xsec_arr[b-1] << " event = " << xsec_i_arr[b-1] << std::endl;
  }

  auto gr = new TGraph( nbin, eg_arr, xsec_i_arr ); gr->SetLineColor(1); gr->SetLineWidth(2);
  gr -> Write();

  return ;

}


