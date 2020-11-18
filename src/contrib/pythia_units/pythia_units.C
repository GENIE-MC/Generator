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


void pythia_units( TString in_file_name  = "gntp.0.ghep.root" ,
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

  TFile out_file ( out_file_name, "RECREATE" ) ;
  out_file.cd() ;

  TTree out_tree( "tau_decay", "#tau decay" ) ;
  
  double length, time ;

  out_tree.Branch( "L", & length, "L/D" );
  out_tree.Branch( "T", & time, "T/D" );


  // Event loop
  for(Long64_t i=0; i < in_tree->GetEntries(); i++) {

    in_tree->GetEntry(i);

    EventRecord & event = *(mcrec->event);

    const Interaction & inter = *( event.Summary() ) ;

    int lep_pdg = inter.FSPrimLeptonPdg() ;

    if ( TMath::Abs( lep_pdg ) == kPdgTau ) { 

      const GHepParticle * tau = event.FinalStatePrimaryLepton() ;

      const GHepParticle * first_daughter = event.Particle( tau -> FirstDaughter() ) ; 

      TLorentzVector distance = * first_daughter -> X4() - * tau -> X4() ;

      length = distance.Vect().Mag() ;
      time = distance.T() ;

      out_tree.Fill() ;

    }
    
    mcrec->Clear() ;
  } // event loop


  out_tree.Write() ;


}
