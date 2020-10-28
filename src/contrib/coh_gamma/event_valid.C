
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
#include "Framework/Conventions/Units.h"

#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"


using namespace genie ;


void event_proc( TTree * in_tree, TFile & ofile ) ;
void event_proc_plot( TTree * in_tree, TFile & ofile ) ;
double LepPhi( const TLorentzVector & probe, TLorentzVector lep ) ;
double GammaPhi( const TLorentzVector & probe, TLorentzVector gamma, TLorentzVector lep ) ;
bool SelEvent( double theta_l, double theta_g, double phi_g ) ;
void GraphIntegratedCOHGammaXsec( double nevts, double xsec, TH1 * h_E_g, TFile & file ) ;

/////////////
// Global constants

double g_width ;
double g_theta_l ;
double g_theta_g ;
double g_phi_g ;
double g_Ev ;
int g_nu ;
int g_target ;

// COH Gamma xsec config set ("Default" or "Consistent")
std::string config = "Default" ;




void event_valid( TString in_file_name  = "../event_gen/numu_e10_12c_500gev.ghep.root" , 
                  TString out_file_name = "",
                  bool loop      = false,
                  int target     = 1000060120,
                  int nu         = 14,
                  double width   = 0.2, 
                  double theta_l = 0.5,
                  double theta_g = 0.5, 
                  double phi_g   = 2.9,
                  double Ev      = 1.0
                ) {

  g_width = width ;
  g_theta_l = theta_l ;
  g_theta_g = theta_g ;
  g_phi_g = phi_g ;
  g_Ev = Ev ;
  g_nu = nu ;
  g_target = target ;

  std::cout << "Starting test with " << g_Ev << " [GeV] probe = " << g_nu << " on target = " << g_target << std::endl;

  // Events < pi/2 and peaked ~0.5
  std::vector<double> theta_angles = {0.3, 0.5, 0.7};
  // Event histogram bowl-shaped so select events on "edges"
  std::vector<double> phi_angle = {-2.9, -2.7, -2.5, -2.3, 2.3, 2.5, 2.7, 2.9};

  if ( out_file_name == "" ) {
    out_file_name = in_file_name ; 
    out_file_name.ReplaceAll( ".ghep.root", 
			      ".evtplots.root" ) ;
  }
    
  TFile in_file( in_file_name ) ; 

  TFile out_file ( out_file_name, "RECREATE" ) ;
  out_file.cd() ;
    
  NtpMCTreeHeader * header = dynamic_cast<NtpMCTreeHeader*> (in_file.Get("header"));

  // Get the GENIE GHEP tree and set its branch address
  TTree * in_tree = dynamic_cast<TTree*> (in_file.Get("gtree"));
  
  if ( loop ) { // Loop over the fixed angles, only xsec plot
    for( auto p : phi_angle ) {
      for( auto t : theta_angles ) {
        g_theta_l = g_theta_g = t; g_phi_g = p ;
        std::cout << "Angle theta_l = theta_g = " << t << " phi_g = " << p << std::endl;
        event_proc( in_tree, out_file ) ;
      }
    }
  } else { // Single angles, make plots in addition to xsec plot
    std::cout << "Angle theta_l = " << g_theta_l << " theta_g = " << g_theta_g << " phi_g = " << g_phi_g << std::endl;
    event_proc_plot( in_tree, out_file ) ; 
  }

  out_file.Close() ;
}


void event_proc_plot( TTree * in_tree, TFile & ofile ) {
  
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
                                                                                                                 
  TH1* h_E_nu = hists["E_nu"] = new TH1D("h_E_nu", "E_{#nu};E_{#nu} [GeV]",
        			       100, 0., 10. ) ; 

  1* h_E_g = hists["E_g"] = new TH1D("h_E_g", "E_{#gamma};E_{#gamma} [GeV]", 20, 0., 1. ) ; 
  h_E_g -> Sumw2( true ) ;

  double nevts = in_tree->GetEntries() ;
  size_t evts = 0;
  size_t sel_evts = 0;
  double xsec ;

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

        xsec = event.XSec()/(1.e-38 * units::cm2) ;
        sel_evts += 1;
	
      } // single gamma

    } // coherent production

      mcrec->Clear() ;
  } // event loop

  GraphIntegratedCOHGammaXsec( nevts, xsec, h_E_g, ofile ) ;

  std::cout << "Processed " << nevts << " and selected " << sel_evts << " events." << std::endl;

  for ( auto & h : hists ) {
    h.second -> Write() ;
    delete h.second ;
  }

}


void event_proc( TTree * in_tree, TFile & ofile ) {

  NtpMCEventRecord * mcrec = 0;
  in_tree->SetBranchAddress("gmcrec", & mcrec);

  TH1* h_E_g = new TH1D("h_E_g", "E_{#gamma};E_{#gamma} [GeV]", 20, 0., 1. ) ;
  h_E_g -> Sumw2( true ) ;

  double nevts = in_tree->GetEntries() ;
  size_t evts = 0;
  size_t sel_evts = 0;
  double xsec ;

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

        if ( !SelEvent( theta_l, theta_g, phi_g ) ) {  mcrec->Clear(); continue; }

         h_E_g -> Fill( gamma.E() ) ;

        xsec = event.XSec()/(1.e-38 * units::cm2) ;
        sel_evts += 1;

      } // single gamma

    } // coherent production

      mcrec->Clear() ;
  } // event loop

  GraphIntegratedCOHGammaXsec( nevts, xsec, h_E_g, ofile ) ;

  std::cout << "Processed " << nevts << " and selected " << sel_evts << " events." << std::endl;

  delete h_E_g ;

  return ;
  
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

  double theta_l_max = g_theta_l + (g_width/2.) ; double theta_l_min = g_theta_l - (g_width/2.) ;
  double theta_g_max = g_theta_g + (g_width/2.) ; double theta_g_min = g_theta_g - (g_width/2.) ;
  double phi_g_max = g_phi_g + (g_width/2.) ; double phi_g_min = g_phi_g - (g_width/2.) ;

  if ( ( theta_l > theta_l_min && theta_l < theta_l_max ) &&
       ( theta_g > theta_g_min && theta_g < theta_g_max ) &&
       ( phi_g   > phi_g_min && phi_g     < phi_g_max ) ) return true;

  return false;
  
}

void GraphIntegratedCOHGammaXsec( double nevts, double xsec, TH1 * h_E_g, TFile & file ) {
  
  // -- request the COH NC Gamma mode 
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("AlvarezRusoSalaCOHGammaPXSec",config) ;
  const Algorithm * algXsec = algf->GetAlgorithm(id);
  const XSecAlgorithmI * fXSec = dynamic_cast<const XSecAlgorithmI *>(algXsec);
   
  Interaction * in = Interaction::COHNC(g_target, g_nu, 22, g_Ev);
  utils::gsl::d4Xsec_dEgdThetaldThetagdPhig func = utils::gsl::d4Xsec_dEgdThetaldThetagdPhig( fXSec, in ) ;

  int nbin = h_E_g -> GetNbinsX() ;
  double xsec_i_arr[nbin];
  double xsec_arr[nbin];
  double eg_arr[nbin];
  double x_err[nbin];
  double y_err[nbin];

  for ( int b = 0; b < nbin; b++ ) {
    double Eg = h_E_g -> GetXaxis() -> GetBinCenter( b ) ;
    double Ni = h_E_g -> GetBinContent( b ) ;
    double wbin = h_E_g -> GetBinWidth( b ) ;
    
    // Calculate xsec from events
    xsec_i_arr[b] =  (1.e3) * xsec * ( Ni / nevts ) / ( wbin * pow( g_width, 3 ) ) ;
    x_err[b] = 0. ; y_err[b] = ( xsec_i_arr[b] / Ni ) * h_E_g -> GetBinError( b ) ; 
    // Get the xsec from the functor
    std::array<double, 4> point = { Eg, g_theta_l, g_theta_g, g_phi_g } ; 
    xsec_arr[b] = (1.e3) * func( point.data() ) ;

    eg_arr[b] = Eg ;
  }

  auto evt_xsec = new TGraphErrors( nbin, eg_arr, xsec_i_arr, x_err, y_err); evt_xsec -> SetLineColor(1); evt_xsec -> SetLineWidth(1); evt_xsec -> SetMarkerStyle(3);
  auto func_xsec = new TGraph( nbin, eg_arr, xsec_arr ); func_xsec -> SetLineColor(2); func_xsec -> SetLineWidth(1); 

  TParticlePDG * tprobe = PDGLibrary::Instance() -> Find( g_nu ) ;
  TParticlePDG * ttgt = PDGLibrary::Instance() -> Find( g_target ) ;

  stringstream description ;
  
  description << g_Ev << " GeV " << tprobe -> GetTitle() << " on " << ttgt -> GetTitle() ;
  description << " #Delta#theta_{l}=" << g_theta_l*TMath::RadToDeg() << "#circ" ;
  description << " #Delta#theta_{#gamma}=" << TMath::RadToDeg()*g_theta_l << "#circ" ;
  description << " #Delta#phi_{#gamma}=" << TMath::RadToDeg()*g_phi_g << "#circ" ;
  description << " ("  << h_E_g -> GetEntries() << "/" << nevts << " evts)" ;
  description << ";E_{#gamma} [GeV];#frac{d^{4}#sigma}{dE_{#gamma}d#theta_{l}d#Omega_{#gamma}} [10^{-41} #frac{cm^{2}}{GeV}]" ;

  std::string evt_title = "Event #sigma: " + description.str() ;
  std::string func_title = "Functor #sigma: " + description.str() ;
  std::string sigma_title = "#sigma: " + description.str() ;

  evt_xsec -> SetTitle( evt_title.c_str() );
  func_xsec -> SetTitle( func_title.c_str() );

  stringstream mg_descrip ;
  mg_descrip << "xsec_angles_" << g_theta_l*TMath::RadToDeg() << "_" << g_theta_g*TMath::RadToDeg() << "_" << g_phi_g*TMath::RadToDeg() ;
  std::string mg_name = mg_descrip.str() ;

  auto mg = new TMultiGraph();
  mg -> Add( func_xsec ) ;
  mg -> Add( evt_xsec ) ;
  mg -> SetTitle( sigma_title.c_str() );
  mg -> Write( mg_name.c_str() ) ;

  return ;

}


