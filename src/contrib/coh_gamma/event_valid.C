
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
double GetXSection( double Eg ) ;


void event_test( TString in_file_name  = "../event_gen/gntp.0.ghep.root" , 
		 TString out_file_name = "" ) {

  
  if ( out_file_name == "" ) {
    out_file_name = in_file_name ; 
    out_file_name.ReplaceAll( ".ghep.root", 
			      ".plots.root" ) ;
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
  TH1* h_sigma = hists["sigma"] = new TH1D("h_sigma", "#sigma;E_{#gamma} [GeV]",
        			       80, 0., 2. ) ; 



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

        evts += 1;

        if ( !SelEvent(lep.Angle( probe.Vect() ), gamma.Angle( probe.Vect() ), GammaPhi( probe, gamma, lep ) ) ) continue;

	h_theta_gamma -> Fill( gamma.Angle( probe.Vect() ) ) ;
	h_theta_lep   -> Fill( lep.Angle( probe.Vect() ) ) ;
	h_theta_recoil   -> Fill( recoil.Angle( probe.Vect() ) ) ;

	h_phi_gamma -> Fill( GammaPhi( probe, gamma, lep ) ) ;			     
	h_phi_lep -> Fill( LepPhi( probe, lep ) ) ;			     


	h_E_l -> Fill( lep.E() ) ;
	h_E_g -> Fill( gamma.E() ) ;
	h_E_nu -> Fill( probe.E() ) ;

        sel_evts += 1;
	
      } // single gamma

    } // coherent production

      mcrec->Clear() ;
  } // event loop

   
  TFile out_file ( out_file_name, "RECREATE" ) ;
  out_file.cd() ;
  
  for ( auto & h : hists ) {
    h.second -> Write() ;
    delete h.second ;
  }

  std::cout << "Processed " << evts << " and selected " << sel_evts << " events." << std::endl;
 
  int steps = 80;
  double xsec_arr[steps];
  double eg_arr[steps];
  double eg = 0;

  for (int i = 0; i < steps; i++  ) { eg += 1; eg_arr[i] = eg/(steps/2); xsec_arr[i] = GetXSection(eg/(steps/2))*(1e3); }

  auto gr = new TGraph( 40, eg_arr, xsec_arr ); gr->SetLineColor(1); gr->SetLineWidth(2);
  gr -> Write();

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

  if ( ( theta_l > 0.0 && theta_l < 0.2 ) &&
       ( theta_g > 0.4 && theta_g < 0.6 ) &&
       ( phi_g > 2.8 && phi_g < 3.0 ) ) return true;

  return false;
  
}

double GetXSection( double Eg ) {
  
  // -- request the COH NC Gamma model
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("AlvarezRusoSalaCOHGammaPXSec","Default");
  const Algorithm * algXsec = algf->GetAlgorithm(id);
  const XSecAlgorithmI * fXSec = dynamic_cast<const XSecAlgorithmI *>(algXsec);
 
  Interaction * in = Interaction::COHNC(1000060120, 14, 22, 2);

  ROOT::Math::IBaseFunctionMultiDim * func = new utils::gsl::d4Xsec_dEgdThetaldThetagdPhig( fXSec, in );

  std::array<double, 4> point = {Eg, 0.1, 0.5, 2.9} ; // { Eg, theta_l, theta_g, phi_g }

  double xsec = (*func)( point.data() ) ;

  if ( xsec < 0 ) std::cout << "Bad Xsection! " << xsec << std::endl;

  //std::cout << "Eg = " << Eg << "  XSec = " << xsec << std::endl;

  return xsec ;

}


