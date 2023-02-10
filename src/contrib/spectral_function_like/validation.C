#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"


using namespace genie ;


void spectral_function_spectra( TString in_file_name  = "gntp.0.ghep.root" ,
                                TString out_file_name = "" ) {


  if ( out_file_name == "" ) {
    out_file_name = in_file_name ;
    out_file_name.ReplaceAll( ".ghep.root",
                              ".plots.root" ) ;
  }

  LOG("Plots", pINFO) << "Reading file: " << in_file_name;
  TFile in_file( in_file_name ) ;

  NtpMCTreeHeader * header = dynamic_cast<NtpMCTreeHeader*> (in_file.Get("header"));

  // Get the GENIE GHEP tree and set its branch address
  TTree * in_tree = dynamic_cast<TTree*> (in_file.Get("gtree"));
  NtpMCEventRecord * mcrec = 0;
  in_tree->SetBranchAddress("gmcrec", & mcrec);


  std::map<std::string, TH1*> hists ;

  // plots we need
  auto temp = hists["Eb_vs_p"] = new TH2D("h_Eb_vs_p", "Spectral function;p_{miss} [GeV];E_{miss} [GeV]",
                                          100, 0., 1., 100, 0., 0.1 ) ;

                                        
  TH2* h_Eb_vs_p = dynamic_cast<TH2*>(temp);
  

  // Event loop
  for(Long64_t i=0; i < in_tree->GetEntries(); i++) {

    in_tree->GetEntry(i);

    EventRecord & event = *(mcrec->event);

    const Interaction & inter = *( event.Summary() ) ;

    const ProcessInfo & proc_info = inter.ProcInfo() ;

    const auto & target = inter.InitState().Tgt() ;

    if ( target.HitNucIsSet() ) {

      if ( ! proc_info.IsMEC() ) {
	   
	const auto * hit_nucleon = event.HitNucleon();
	
	h_Eb_vs_p->Fill( hit_nucleon->P4()->P(), hit_nucleon->RemovalEnergy() );
	
	if ( hit_nucleon->RemovalEnergy() < 0.005 ) 
	  cout << event << endl;
      }
    }
    
    mcrec->Clear() ;
  } // event loop


  LOG("Plots", pINFO) << "Writing to file: " << out_file_name;
  TFile out_file ( out_file_name, "RECREATE" ) ;
  out_file.cd() ;

  for ( auto & h : hists ) {
    SLOG("Plots", pINFO) << "Plot: " << h.second->GetTitle();
    h.second -> Write() ;
    delete h.second ;
  }
 
}
