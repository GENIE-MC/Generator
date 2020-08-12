
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
#include "Framework/ParticleData/PDGUtils.h" 

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
  TH1* h_E_N = hists["E_N"] = new TH1D("h_E_N", "E_{N};E_{N} [GeV]",
				       80, 0., 6. ) ; 

  TH1* h_theta_N = hists["theta_N"] = new TH1D("h_theta_n", "#theta_{N};#theta_{N} [rad]",
					       80, 0., TMath::Pi() ) ; 


  // Event loop
  for(Long64_t i=0; i < in_tree->GetEntries(); i++) {
 
    in_tree->GetEntry(i);
    
    EventRecord & event = *(mcrec->event);
    
    const Interaction & inter = *( event.Summary() ) ;
    
    const ProcessInfo & proc_info = inter.ProcInfo() ;
    
    if ( proc_info.IsCoherentElastic() ) {
    
      if ( proc_info.IsDarkNeutralCurrent() ) {
      
	const GHepParticle & N = * event.Particle(2) ; 
	
	const GHepParticle * final_neutrino = nullptr ; 
	const GHepParticle * mediator = nullptr ; 
	std::vector<const GHepParticle*> final_products ;
	
	const GHepParticle * temp = event.Particle( N.FirstDaughter() ) ;
	
	// setting the final products might not be easy for ever
	// let's try to get it right once and for all
	if ( N.LastDaughter() - N.FirstDaughter() == 1 ) { 
	  // the dark neutrino decays only in two bodies
	  // so one is the neutrino
	  // the other is the dark mediator
	  if ( pdg::IsNeutrino( TMath::Abs( temp -> Pdg() ) ) ) { 
	    final_neutrino = temp ;
	    mediator = event.Particle( N.LastDaughter() ) ;
	  }
	  else { 
	    mediator = temp ;
	    final_neutrino = event.Particle( N.LastDaughter() ) ;
	  }
	  
	  // the final products are then the daughters of the mediator
	  for ( unsigned int i = mediator -> FirstDaughter() ;
		i <= mediator -> LastDaughter() ; ++i ) {
	    
	    final_products.push_back( event.Particle( i ) ) ;
	  }

	}	  
	else {
	  // otherwise there is a neutrino and there is no mediator
	  // so the neutrino is among the decay products fo the dark neutrno, 
	  // all the rest is products
	  // impossible to say which one is the main neutrino, just pick the first in that case
	  
	  for ( unsigned int i = N.FirstDaughter() ;
		i <= N.LastDaughter() ; ++i ) {
	    
	    temp = event.Particle( i ) ;
	    if ( ! final_neutrino ) {
	      if ( pdg::IsNeutrino( TMath::Abs( temp -> Pdg() ) ) ) {
		final_neutrino = temp ;
		continue ;
	      }
	    }
	    final_products.push_back( temp ) ;
	  }

	}  // the neutrino decays in other than 2 daughters

	// now we have all the particles identified and we can fill the hists
	// note that the mediator pointer might be 0 as it does not necessarily exist
	
	const TLorentzVector & probe = * event.Probe() -> P4() ;
	
	const TLorentzVector & p4_N =  * N.P4()  ; 
	const TLorentzVector & p4_recoil = * event.Particle(3) ->P4() ; 

	
	h_E_N -> Fill( p4_N.E() ) ;
	h_theta_N -> Fill( p4_N.Angle( probe.Vect() ) ) ;
	
      } // dark neutral current 

    } // coherent elastic

      mcrec->Clear() ;
  } // event loop

   
  TFile out_file ( out_file_name, "RECREATE" ) ;
  out_file.cd() ;
  
  for ( auto & h : hists ) {
    h.second -> Write() ;
    delete h.second ;
  }

}


