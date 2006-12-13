//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Prints the events in the input GENIE/ROOT file
//
// Usage:
//   shell% genie
//   genie[0] .x printEvents.C("/path/to/my/genie_file.root");
//
// Costas Andreopoulos, Dec 13, 2006
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void printEvents(const char * filename)
{
 //-- open the ROOT file and get the event TTree 
 TFile inpfile(filename,"READ");
 TTree * tree = dynamic_cast <TTree *> (inpfile.Get("gtree"));

 //-- set the branch address
 genie::NtpMCEventRecord * mcrec = 0;
 tree->SetBranchAddress("gmcrec", &mcrec);

 //-- loop over TTree NtpMC records
 for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);

    genie::NtpMCRecHeader rec_header = mcrec->hdr;
    genie::EventRecord &  event      = *(mcrec->event);

    cout << rec_header;
    cout << event;

    mcrec->Clear();
 }
 inpfile.Close();
}

