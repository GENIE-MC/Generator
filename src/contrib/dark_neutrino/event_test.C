
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

  TH1* h_T_T = hists["T_T"] = new TH1D("h_T_T", "T_{T};T_{T} [GeV]",
				       100, 0., 0.001 ) ; 
  
  TH1* h_E_vis = hists["E_vis"] = new TH1D("h_E_vis", "E_{vis};E_{vis} [GeV]",
					   100, 0., 6. ) ; 

  TH1* h_vis_dec = hists["vis_dec"] = new TH1D("h_vis_dec", "Visible Decay;visible",
					       2, -0.5, 1.5 ) ; 


  TH1* h_E_prod = hists["E_prod"] = new TH1D("h_E_prod", "E_{prod};E_{prod} [GeV]",
					     100, 0., 1. ) ; 


  TH1* h_theta_N = hists["theta_N"] = new TH1D("h_theta_n", "#theta_{N};#theta_{N} [rad]",
					       80, 0., TMath::Pi() ) ; 

  TH1* h_N_prod = hists["N_prod"] = new TH1D("h_N_prod", "N_{Prod};N_{Prod}",
					     10, -0.5, 9.5 ) ; 


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
	
	h_N_prod -> Fill( final_products.size() ) ;

	const TLorentzVector & probe = * event.Probe() -> P4() ;
	
	const TLorentzVector & p4_N =  * N.P4()  ; 
	const TLorentzVector & p4_recoil = * event.Particle(3) ->P4() ; 

	
	h_E_N -> Fill( p4_N.E() ) ;
	
	double t_t = p4_recoil.E() - p4_recoil.Mag() ;
	h_T_T -> Fill( t_t ) ;

	h_theta_N -> Fill( p4_N.Angle( probe.Vect() ) ) ;
	
	

	double vis_e = t_t ; 
	double prod_e = 0. ;
	bool vis_decay = false ;
	for ( const auto & p : final_products ) {
	  
	  if ( p -> Pdg() == 22 )  vis_e += p -> P4() -> E() ;
	  else if ( p -> Charge() != 0. ) { 
	    vis_decay = true ;
	    vis_e += p -> P4() -> E() ;
	  }
	  
	  prod_e += p -> P4() -> E() ;
	
	} // finale products loop

	if ( mediator ) {
	  // evaluate the lenght of the first and second day
	  
	  TLorentzVector Delta = (* mediator -> X4()) - (* N.X4()) ;
	  double total_distance = Delta.Vect().Mag() ;
	  double delay = Delta.T() ;
	  
	  std::cout << "Dark neutrino: " << total_distance << " fm in " << delay << " 10^-24 s" << std::endl ;

	  Delta = (*final_products[0] -> X4()) - (* mediator -> X4()) ;
	  total_distance = Delta.Vect().Mag() ;
	  delay = Delta.T() ;
	  
	  std::cout << "mediator: " << total_distance << " fm in " << delay << " 10^-24 s" << std::endl ;

	}
	
	TLorentzVector Delta = (*final_products[0] -> X4()) - (* N.X4() ) ;
	double total_distance = Delta.Vect().Mag() ;
	double delay = Delta.T() ;
	
	std::cout << "Total: " << total_distance << " fm in " << delay << " 10^-24 s" << std::endl ;
	
	h_vis_dec -> Fill( vis_decay ? 1 : 0 ) ;
	h_E_vis -> Fill( vis_e ) ;
	h_E_prod -> Fill( prod_e ) ;
	
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


