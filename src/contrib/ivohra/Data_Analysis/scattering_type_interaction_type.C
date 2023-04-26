//____________________________________________________________________________
/*
\brief    Plots a 2D histogram depicting the distribution of the scattering
          type and the interaction type of each event.
          
          To see the documentation for this example, visit: 
          https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=349
          
\author   Ishaan Vohra <ivohra@exeter.edu / ishaanklv@gmail.com>
          Phillips Exeter Academy
\created  August 16, 2022
*/
//____________________________________________________________________________

#include <iostream>
#include <string>
#include "GHepParticle.h"
#include "ScatteringType.h"
#include "ProcessInfo.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TAxis.h"
#include "NtpMCTreeHeader.h"
#include "NtpMCEventRecord.h"
#include "EventRecord.h"
#include "BaryonResUtils.h"
#include "GHepParticle.h"
#include "PDGUtils.h"
#include "GHepStatus.h"
#include "PDGCodes.h"

using namespace genie;

void scattering_type_interaction_type(string filename = "/hepstore/ivohra/miniboone1.ghep.root")
{

  // Open the GHEP ROOT file

  TFile infile(filename.c_str());

  // Get the tree header

  NtpMCTreeHeader *header =
      dynamic_cast<NtpMCTreeHeader *>(infile.Get("header"));

  // Get the GENIE GHEP tree and set its branch address

  TTree *tree = dynamic_cast<TTree *>(infile.Get("gtree"));
  NtpMCEventRecord *mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  // Make the histogram

  auto myHist = new TH2D("h1", "Distribution of Interaction Type and Scattering Type; ", 20, -0.5, 19.5, 3, 1.5, 4.5);

  // Obtain the x-axis bin centers and label each bin

  for (int i = 1; i <= myHist->GetNbinsX(); i++)
  {
    auto binCenter = myHist->GetXaxis()->GetBinCenter(i);

    string s = ScatteringType::AsString(ScatteringType_t(binCenter));
    myHist->GetXaxis()->SetBinLabel(i, s.c_str());
  }

  // Label the y-axis bins

  myHist->GetYaxis()->SetBinLabel(1, "CC");
  myHist->GetYaxis()->SetBinLabel(2, "NC");
  myHist->GetYaxis()->SetBinLabel(3, "Mix");

  // Event loop

  for (Long64_t i = 0; i < tree->GetEntries(); i++)
  {
    tree->GetEntry(i);

    EventRecord &event = *(mcrec->event);
    const Interaction *in = event.Summary();
    const ProcessInfo &proc = in->ProcInfo();

    // Fill the bins

    myHist->Fill(proc.ScatteringTypeId(), proc.InteractionTypeId());
    mcrec->Clear();
  }

  // Normalization and stat box removal

  TH2 *normh = (TH2D *)(myHist->Clone("normh"));
  normh->Scale(1. / normh->GetEntries());
  normh->SetStats(false);

  // Create the output file

  TFile *outfile = new TFile("scattering_type_interaction_type.root", "RECREATE");
  outfile->cd();
  normh->Write();
  outfile->Close();
}
