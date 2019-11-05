#include "TString.h" 
#include "TFile.h"

#include "Framework/Ntuple/NtpMCTreeHeader.h" 


int make_plots( TString file_name ) {

  
  TFile in_file( file_name );

  NtpMCTreeHeader * header = dynamic_cast<NtpMCTreeHeader*>( file.Get("header") );

  // Get the GENIE GHEP tree and set its branch address
  TTree * tree = dynamic_cast<TTree*> ( file.Get("gtree") );
  
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress( "gmcrec", &mcrec);

  // Event loop
  for(Long64_t i=0; i < tree->GetEntries(); i++){
    
    tree->GetEntry(i);
    // print-out the event
    
    EventRecord & event = *(mcrec->event);
    
    


    mcrec->Clear();
  }

  return tree -> GetEntries() ;

}
