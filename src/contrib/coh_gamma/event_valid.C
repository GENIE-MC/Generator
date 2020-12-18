/////////////////////////////////////////////////////////
//
// To run use `genie -l event_valid.C`
//
//  or `genie -l "event_valid.C( int nbins, TString input_file, TString output_file )"`
//   - e.g: `genie -l "event_valid.C( 5, \"vmu_12c_1gev.root\", \"test.root\" )"`
//
// Parameters:
//
//  Input file:  in_file_name (line 76)
//  Output file: out_file_name (line 77)
//
//  Bin number: nbins (line 75)
//    - bins for theta_l = theta_g = nbins
//    - bins for phi_g = nbins * 2
//    - bins for E_g = nbins * 6
//
//  J. Sensenig 12/2020
//
//////////////////////////////////////////////////////

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


void SetValues( TTree * in_tree ) ;
void event_proc( TTree * in_tree, TFile & ofile ) ;
double LepPhi( const TLorentzVector & probe, TLorentzVector lep ) ;
double GammaPhi( const TLorentzVector & probe, TLorentzVector gamma, TLorentzVector lep ) ;
void GetLimits( std::pair<double,double> &theta_l ) ;
void ComputeCOHGammXSec( double xsec, TH3 * h_Eg_thetag_phig );
void GraphIntegratedCOHGammaXsec( int nbin, double nevts, double *xsec_i_arr, double *xsec_arr, double *eg_arr, double *x_err, double *y_err );

/////////////////////
// Global constants

int g_nbins ;
double g_theta_l_width ;
double g_theta_l ;
double g_theta_g ;
double g_theta_g_max ;
double g_phi_g ;
double g_phi_g_max ;
double g_phi_g_min ;
double g_Ev ;
double g_Eg_max ;
double g_nevts;
int g_nu ;
int g_target ;

// COH Gamma xsec config set ("Default" or "Consistent")
std::string config = "Default" ;


void event_valid( int nbins  = 5, //use 10 bins
                  TString in_file_name  = "../12c/vmu/0.5/vmu_on_1000060120_COHGAMMA_0.5GeV.gntp.root" , 
                  TString out_file_name = "test_line_fast_out_05gev.root"
                ) {

  g_nbins = nbins ;

  if ( out_file_name == "" ) {
    out_file_name = in_file_name ; 
    out_file_name.ReplaceAll( ".ghep.root", ".evtplots.root" ) ;
  }
    
  TFile in_file( in_file_name ) ; 

  if ( !in_file.IsOpen() ) {

  std::cout << "File " << in_file_name << " not open! Exiting!" << std::endl ;
  return ;

  }

  TFile out_file ( out_file_name, "RECREATE" ) ;
  out_file.cd() ;
    
  // Get the GENIE GHEP tree and set its branch address
  TTree * in_tree = dynamic_cast<TTree*> (in_file.Get("gtree"));

  std::cout << "Extracting and setting values" << std::endl;

  // Get the data ranges
  SetValues( in_tree ) ;
  
  int loop_cnt = 0;

  for ( int  i = 0; i < g_nbins; i++ ) { // theta_l

    loop_cnt += 1;
    std::cout << "Loop " << loop_cnt << "/" << g_nbins << " with theta_l = " << g_theta_l << std::endl;

    event_proc( in_tree, out_file ) ;
    g_theta_l += g_theta_l_width ;

  }

  out_file.Close() ;
}

void SetValues( TTree * in_tree ) {

  // Grab interaction details from first event
  NtpMCEventRecord * mcrec = 0;
  in_tree->SetBranchAddress("gmcrec", & mcrec);
  in_tree->GetEntry(0);
  EventRecord & event = *(mcrec->event);

  g_nu = event.Probe() -> Pdg() ;
  g_Ev = event.Probe() -> E() ;
  g_target = event.Particle(4) -> Pdg() ;

  mcrec->Clear();

  double Eg_max      = 0. ;
  double theta_l_max = 0. ;
  double theta_g_max = 0. ;
  double phi_g_max   = 0. ;
  double phi_g_min   = 0. ;

  // Get the max/min limits from the data
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

        double Eg = gamma.E() ;
        double theta_l = lep.Angle( probe.Vect() ) ;
        double theta_g = gamma.Angle( probe.Vect() ) ;
        double phi_g = GammaPhi( probe, gamma, lep ) ;

        if ( Eg > Eg_max )           Eg_max = Eg ;
        if ( theta_l > theta_l_max ) theta_l_max = theta_l ;
        if ( theta_g > theta_g_max ) theta_g_max = theta_g ;
        if ( phi_g > phi_g_max )     phi_g_max = phi_g ;
        if ( phi_g < phi_g_min )     phi_g_min = phi_g ;

      }
    }
    mcrec->Clear();
  }
 
  // Long tail avoid for now
  if ( theta_l_max > TMath::PiOver2() ) theta_l_max = TMath::PiOver2() ;  
  if ( theta_g_max > TMath::PiOver2() ) theta_g_max = TMath::PiOver2() ;  

  std::cout << " Ev = "     << g_Ev 
            << " Probe = "  << g_nu 
            << " Target = " << g_target << std::endl;

  std::cout << " Eg max = "      << Eg_max 
            << " theta_l max = " << theta_l_max 
            << " theta_g max = " << theta_g_max 
            << " phi_g max = "   << phi_g_max
            << " phi_g min = "   << phi_g_min << std::endl;
  
  g_theta_l_width = theta_l_max / g_nbins ;
  g_theta_l       = g_theta_l_width ;
  g_theta_g_max   = theta_g_max;
  g_phi_g_max     = phi_g_max;
  g_phi_g_min     = phi_g_min;
  g_Eg_max        = Eg_max ;

}



void event_proc( TTree * in_tree, TFile & ofile ) {

  std::pair<double,double> theta_l_lim ; 

  GetLimits( theta_l_lim ) ;

  NtpMCEventRecord * mcrec = 0;
  in_tree->SetBranchAddress("gmcrec", & mcrec);

  TH3* h_Eg_thetag_phig = new TH3D("h_Eg_thetag_phig", "E_{#gamma} #theta_{#gamma} #phi_{#gamma};E_{#gamma} [GeV];#theta_{#gamma} [rad];#phi_{#gamma} [rad]",
                                    (g_nbins*6), 0., g_Eg_max, g_nbins, 0, g_theta_g_max, (g_nbins*2), g_phi_g_min, g_phi_g_max ) ;

  g_nevts = in_tree->GetEntries() ;
  size_t evts = 0;
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

        if ( theta_l < theta_l_lim.first || theta_l > theta_l_lim.second )  {  mcrec->Clear(); continue; }
        
        h_Eg_thetag_phig -> Fill( gamma.E(), theta_g, phi_g ) ;

        xsec = event.XSec()/(1.e-38 * units::cm2) ;

      } // single gamma

    } // coherent production

      mcrec->Clear() ;
  } // event loop

  ComputeCOHGammXSec( xsec, h_Eg_thetag_phig ) ;

  h_Eg_thetag_phig -> Write();

  std::cout << "Processed " << g_nevts << " events." << std::endl;

  delete h_Eg_thetag_phig ;

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


void GetLimits( std::pair<double,double> &theta_l ) {
 
  theta_l = std::make_pair( g_theta_l - (g_theta_l_width/2.), g_theta_l + (g_theta_l_width/2.) ) ;

}


void ComputeCOHGammXSec( double xsec, TH3 * h_Eg_thetag_phig ) {

  // -- request the COH NC Gamma mode
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("AlvarezRusoSalaCOHGammaPXSec",config) ;
  const Algorithm * algXsec = algf->GetAlgorithm(id);
  const XSecAlgorithmI * fXSec = dynamic_cast<const XSecAlgorithmI *>(algXsec);

  Interaction * in = Interaction::COHNC(g_target, g_nu, 22, g_Ev);
  utils::gsl::d4Xsec_dEgdThetaldThetagdPhig func = utils::gsl::d4Xsec_dEgdThetaldThetagdPhig( fXSec, in ) ;

  int nbins_Eg     = h_Eg_thetag_phig -> GetNbinsX(); // X = E_g
  int nbins_thetag = h_Eg_thetag_phig -> GetNbinsY(); // Y = theta_g
  int nbins_phig   = h_Eg_thetag_phig -> GetNbinsZ(); // Z = phi_g

  std::cout << "Processing Eg bins " << nbins_Eg << " theta_g bins " <<  nbins_thetag << " phi_g bins " <<  nbins_phig << std::endl;

  // Loop over bins on each axis ignoring the under/over flow bins
  for ( int b_theta = 1; b_theta < (nbins_thetag+1); b_theta++ ) { // theta_g
    for ( int b_phi = 1; b_phi < (nbins_phig+1); b_phi++ ) {     // phi_g

      double theta_g       = h_Eg_thetag_phig -> GetYaxis() -> GetBinCenter( b_theta );
      double theta_g_width = h_Eg_thetag_phig -> GetYaxis() -> GetBinWidth( b_theta );
      double phi_g         = h_Eg_thetag_phig -> GetZaxis() -> GetBinCenter( b_phi );
      double phi_g_width   = h_Eg_thetag_phig -> GetZaxis() -> GetBinWidth( b_phi );

      double xsec_i_arr[nbins_Eg];
      double xsec_arr[nbins_Eg];
      double eg_arr[nbins_Eg];
      double x_err[nbins_Eg];
      double y_err[nbins_Eg];

      // Project out Eg so we can get the bin error
      TH1D * h_E_g = h_Eg_thetag_phig -> ProjectionX( "h_E_g", b_theta, b_theta, b_phi, b_phi, "e" );

      double cnt = 0;

      for ( int b = 1; b < (nbins_Eg+1); b++ ) { // E_g

        double Eg       = h_Eg_thetag_phig -> GetXaxis() -> GetBinCenter( b );
        double Eg_width = h_Eg_thetag_phig -> GetXaxis() -> GetBinWidth( b );
        double Ni       = h_Eg_thetag_phig -> GetBinContent( b, b_theta, b_phi );

        // Calculate xsec from events
        xsec_i_arr[b-1] = ( 1.e3 ) * xsec * ( Ni / g_nevts ) / ( Eg_width * g_theta_l_width * theta_g_width * phi_g_width );

        x_err[b] = Eg_width / 2.;
        y_err[b] = ( xsec_i_arr[b-1] / Ni ) * h_E_g -> GetBinError( b );

        // Get the xsec from the functor at given 4d point
        std::array<double, 4> point = {Eg, g_theta_l, theta_g, phi_g};
        xsec_arr[b-1] = ( 1.e3 ) * func( point.data() ) ;

        eg_arr[b-1] = Eg;

        cnt += Ni;

      } // E_g

      if ( cnt < 1 ) continue; // Skip plotting if no events in bin

      g_theta_g = theta_g; g_phi_g = phi_g;
      GraphIntegratedCOHGammaXsec( nbins_Eg, cnt, xsec_i_arr, xsec_arr, eg_arr, x_err, y_err );

    } // phi_g
  }   // theta_g

}


void GraphIntegratedCOHGammaXsec( int nbin, double nevts, double *xsec_i_arr, double *xsec_arr, double *eg_arr, double *x_err, double *y_err ) {
  
  auto evt_xsec = new TGraphErrors( nbin, eg_arr, xsec_i_arr, x_err, y_err); 
  evt_xsec -> SetLineColor(1); evt_xsec -> SetLineWidth(1); evt_xsec -> SetMarkerStyle(8); evt_xsec -> SetMarkerSize(0.5);

  auto func_xsec = new TGraph( nbin, eg_arr, xsec_arr ); 
  func_xsec -> SetLineColor(2); func_xsec -> SetLineWidth(1); func_xsec -> SetMarkerStyle(8); func_xsec -> SetMarkerSize(0.3);

  TParticlePDG * tprobe = PDGLibrary::Instance() -> Find( g_nu ) ;
  TParticlePDG * ttgt = PDGLibrary::Instance() -> Find( g_target ) ;

  stringstream description ;
  
  description << g_Ev << " GeV #"    << tprobe -> GetTitle() << " on " << ttgt -> GetTitle() ;
  description << " #theta_{l}="      << (int)( TMath::RadToDeg()*g_theta_l ) << "#circ" ;
  description << " #theta_{#gamma}=" << (int)( TMath::RadToDeg()*g_theta_g ) << "#circ" ;
  description << " #phi_{#gamma}="   << (int)( TMath::RadToDeg()*g_phi_g )     << "#circ" ;
  description << " ("  << nevts << "/" << g_nevts << " evts)" ;
  description << ";E_{#gamma} [GeV];#frac{d^{4}#sigma}{dE_{#gamma}d#theta_{l}d#Omega_{#gamma}} [10^{-41} #frac{cm^{2}}{GeV}]" ;

  std::string evt_title = "Event #sigma: " + description.str() ;
  std::string func_title = "Functor #sigma: " + description.str() ;
  std::string sigma_title = "#sigma: " + description.str() ;

  evt_xsec -> SetTitle( evt_title.c_str() );
  func_xsec -> SetTitle( func_title.c_str() );

  stringstream mg_descrip ;
  mg_descrip << "xsec_angles_" << (int)( g_theta_l*TMath::RadToDeg() ) 
                        << "_" << (int)( g_theta_g*TMath::RadToDeg() ) 
                        << "_" << (int)( g_phi_g*TMath::RadToDeg() ) ;
  std::string mg_name = mg_descrip.str() ;

  auto mg = new TMultiGraph();
  mg -> Add( func_xsec ) ;
  mg -> Add( evt_xsec ) ;
  mg -> SetTitle( sigma_title.c_str() );
  mg -> Write( mg_name.c_str() ) ;

  return ;

}


