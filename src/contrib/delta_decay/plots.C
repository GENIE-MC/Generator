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


int FindDelta( const EventRecord & event ) ;
int IsSinglePion( const EventRecord & event ) ; 




int make_plots( TString file_name, TString flag, TString dir ) {

  
  TFile in_file( dir + '/' + file_name );

  //NtpMCTreeHeader * header = dynamic_cast<NtpMCTreeHeader*>( file.Get("header") );

  // Get the GENIE GHEP tree and set its branch address
  TTree * tree = dynamic_cast<TTree*> ( in_file.Get("gtree") );

  TH1D h_ct( "h_cos_theta", flag + " cos(#theta) distribution;Cos(#theta)", 100, -1., 1. ) ;
  TH1D h_ct_lab( "h_cos_theta_lab", flag + " cos(#theta) distribution;Cos(#theta)", 100, -1., 1. ) ;

  TH1D h_theta( "h_theta", flag + " #theta distribution;#theta [rad]", 100, 0., TMath::Pi() ) ;
  TH1D h_phi  ( "h_phi",   flag + " #varphi distribution;#varphi [rad]", 100, 0., TMath::Pi() ) ;
  
  TH2D h_tp( "h_theta_phi", flag + " (#theta, #varphi) distribution;#theta [rad];#varphi [rad]", 
	     50, 0., TMath::Pi(),
	     50, 0., TMath::Pi() ) ;

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress( "gmcrec", &mcrec);

  // Event loop
  for(Long64_t i=0; i < tree->GetEntries(); i++){
    
    tree->GetEntry(i);
    // print-out the event
    
    EventRecord & event = *( mcrec->event );
    
    const Interaction & inter = *( event.Summary() ) ;

    const ProcessInfo & proc_info = inter.ProcInfo() ;

    if ( proc_info.IsResonant() ) { 

      const XclsTag & xcl_tag = inter.ExclTag() ;

      Resonance_t res = xcl_tag.Resonance() ;

      if ( utils::res::IsDelta( res )  ) {

	int delta_pos = FindDelta( event ) ;

	if ( delta_pos < 0 ) { 
	  cerr << "found delta event with no delta" << std::endl ;
	  cerr << "Event " << i << std::endl ;
	  mcrec->Clear();
	  continue  ; 
	}

	GHepParticle * delta = event.Particle( delta_pos ) ; 

	uint n_d = delta -> LastDaughter() - delta -> FirstDaughter() + 1 ; 

	if ( n_d == 2 ) {

	  GHepParticle * d1 = event.Particle( delta -> FirstDaughter() ) ;
	  GHepParticle * d2 = event.Particle( delta -> LastDaughter() ) ;

	  if ( ( pdg::IsPion( d1 -> Pdg() )    || pdg::IsPion( d2 -> Pdg() ) ) && 
	       ( pdg::IsNucleon( d1 -> Pdg() ) || pdg::IsNucleon( d2 -> Pdg() ) ) ) {
	    
	    GHepParticle * pion = pdg::IsPion( d1 -> Pdg() ) ? d1 : d2 ;
	    //GHepParticle * nucl = pdg::IsNucleon( d1 -> Pdg() ) ? d1 : d2 ;
	    
	    const TLorentzVector & delta_mom = *( delta -> P4() ) ;
	    
	    TLorentzVector pion_mom = * ( pion -> P4() ) ;
	    //TLorentzVector nucl_mom = * ( nucl -> P4() ) ;
	    
	    TLorentzVector in_lep_p4( * (event.Probe()-> P4() ) ) ;
	    
	    TLorentzVector out_lep_p4 = *( event.FinalStatePrimaryLepton() -> P4() ) ;
	    
	    TLorentzVector q = in_lep_p4 - out_lep_p4 ;
	    
	    double theta_lab = pion_mom.Angle( q.Vect() ) ;
	    h_ct_lab.Fill( TMath::Cos( theta_lab ) ) ;

	    // move in the Delta CM

	    in_lep_p4.Boost( - delta_mom.BoostVector() ) ;
	    out_lep_p4.Boost( - delta_mom.BoostVector() ) ;
	    
	    q = in_lep_p4 - out_lep_p4 ;
	    
	    pion_mom.Boost( - delta_mom.BoostVector() );  // this gives us the pion in the Delta reference frame
	    

	    TVector3 pion_dir = pion_mom.Vect().Unit() ;
	    TVector3 z_axis = q.Vect().Unit() ;

	    TVector3 y_axis = in_lep_p4.Vect().Cross( out_lep_p4.Vect() ).Unit() ;
	    TVector3 x_axis = y_axis.Cross(z_axis);

	    TVector3 pion_perp( z_axis.Cross( pion_dir.Cross( z_axis ).Unit() ) ) ;
    
	    double phi = pion_perp.Angle( x_axis ) ;  // phi
	    
	    double theta = pion_dir.Angle( z_axis ) ;
	    
	    h_ct.Fill( TMath::Cos(theta) ) ;
	    h_theta.Fill( theta ) ;
	    h_phi.Fill( phi ) ;
	    h_tp.Fill( theta, phi ) ;

	  }  // if there is a pion and a nucleon

	} // if the delta has 2 daughters


	
	
      } // is Delta 

      
      
    } // is resonant process


    // int pion_pos = IsSinglePion( event ) ;
    // if ( pion_pos > 0 ) {
      
    //   TLorentzVector in_lep_p4( * (event.Probe()-> P4() ) ) ;
      
    //   TLorentzVector out_lep_p4 = *( event.FinalStatePrimaryLepton() -> P4() ) ;
      
    //   TLorentzVector q = in_lep_p4 - out_lep_p4 ;
      
    //   TLorentzVector pion_mom( * event.Particle( pion_pos ) -> P4() ) ;

    //   double theta = pion_mom.Angle( q.Vect() ) ;
    //   h_ct_lab.Fill( TMath::Cos( theta ) ) ;
		     
    // }  //  we have a single pion event 

      
    mcrec->Clear();
  }

  TFile out_file( "result_" + file_name, "RECREATE" ) ;
  out_file.cd() ;
  
  h_ct.Write() ;
  h_ct_lab.Write() ;
  h_theta.Write() ;
  h_phi.Write() ;
  h_tp.Write() ;

  return tree -> GetEntries() ;

}


int FindDelta( const EventRecord & event ) {

  TObjArrayIter iter(&event);
  GHepParticle * p = 0;
  GHepParticle * d = 0;
  
  int p_id = -1 ;
  // loop over event particles
  while( (p = dynamic_cast<GHepParticle *>(iter.Next())) ) {
    
    p_id ++ ; 
  
    int pdgc = p->Pdg();
    int status = p->Status();

    if( status != kIStDecayedState ) continue;

    if ( ! utils::res::IsDelta( utils::res::FromPdgCode(pdgc) ) ) continue ;

    bool has_nucleon = false ; 
    for ( int i = p -> FirstDaughter() ;
	  i <= p -> LastDaughter() ; ++i ) {
      
      d = event.Particle( i ) ;
      if ( pdg::IsNucleon( d -> Pdg() ) ) {
	has_nucleon = true ;
	break ;
      }
      
    }

    if ( has_nucleon ) return p_id ;
     
  }

  return -1 ;
}

int IsSinglePion( const EventRecord & event ) {

  TObjArrayIter iter(&event);
  GHepParticle * p = 0;
    
  unsigned int n_pions = 0 ;
  unsigned int n_protons = 0 ;
  unsigned int n_leptons = 0 ;
  unsigned int n_others = 0 ;

  int pion_pos = -1 ; 
  int prot_pos = -1 ;
  int pos = -1 ;
  // loop over event particles
  while( (p = dynamic_cast<GHepParticle *>(iter.Next())) ) {
    
    pos ++ ; 

    int pdgc = p->Pdg();
    int status = p->Status();

    if( status != kIStDecayedState ) continue;

    if      ( genie::pdg::IsChargedLepton( pdgc ) ) n_leptons ++  ;
    else if ( genie::pdg::IsNucleon( pdgc ) ) {
      if ( genie::pdg::IsProton( pdgc ) ) {
	if ( p -> KinE() > 0.1 ) {
	  prot_pos = pos ;
	  n_protons ++ ;
	}
	
      }
    }

    else if ( genie::pdg::IsPion( pdgc ) ) { 
      if ( TMath::Abs(pdgc) == genie::kPdgPiP ) {
	pion_pos = pos ;
	n_pions++ ;
      }
    }
    else if ( genie::pdg::IsParticle( pdgc ) ) {
      if ( p -> Charge() != 0 ) 
	n_others ++ ;
    }
    
    
  }  // particle loop
  
  if ( n_leptons > 1 ) return -1 ;
  if ( n_others > 0 ) return -1 ;
  if ( n_protons != 1 ) return -1 ;
  if ( n_pions != 1 ) return -1 ;
  
  TLorentzVector delta_mom = (* event.Particle( prot_pos ) -> P4()) + (* event.Particle( pion_pos ) -> P4()) ;
  double w = delta_mom.M() ; 

  if ( w < 1.4 ) return pion_pos ;
  return -1 ; 
}
