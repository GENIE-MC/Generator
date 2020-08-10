
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

using namespace genie ;





void event_test( TString in_file_name  = "gntp.0.ghep.root" , 
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

	h_theta_gamma -> Fill( lep.Angle( probe.Vect() ) ) ;

	
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

}

