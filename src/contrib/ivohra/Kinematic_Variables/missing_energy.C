//____________________________________________________________________________
/*
\brief    Assuming all neutral particles are unobserved, plots a 1D histogram
          depicting the distribution of the total final state missing energy
          of each event.
          
          To see the documentation for this example, visit: 
          https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=349
          
\author   Ishaan Vohra <ivohra@exeter.edu / ishaanklv@gmail.com>
          Phillips Exeter Academy
\created  August 16, 2022
*/
//____________________________________________________________________________

#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1D.h"
#include "TH2.h"
#include "GHepParticle.h"
#include "ScatteringType.h"
#include "ProcessInfo.h"
#include "TObjArray.h"
#include "PDGLibrary.h"
#include "TString.h"
#include "NtpMCTreeHeader.h"
#include "NtpMCEventRecord.h"
#include "EventRecord.h"
#include "BaryonResUtils.h"
#include "PDGUtils.h"
#include "GHepStatus.h"
#include "PDGCodes.h"
#include "TDatabasePDG.h"

using namespace genie;

void missing_energy(string filename = "/hepstore/ivohra/miniboone1.ghep.root")
{

  // Open the GHEP/ROOT file

  TFile infile(filename.c_str());

  // Get the tree header

  NtpMCTreeHeader *header =
      dynamic_cast<NtpMCTreeHeader *>(infile.Get("header"));

  // Get the GENIE GHEP tree and set its branch address

  TTree *tree = dynamic_cast<TTree *>(infile.Get("gtree"));
  NtpMCEventRecord *mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  // Make the histogram

  auto myHist = new TH1D("h1", "Missing Energy (Neutral Particle Energy);Energy (GeV);Probability", 100, 0., 5.);

  // Event loop

  for (Long64_t i = 0; i < (tree->GetEntries()); i++)
  {
    tree->GetEntry(i);

    EventRecord &event = *(mcrec->event);
    TObjArrayIter iter(&event);
    GHepParticle *p = 0;

    double e0Sum = 0;

    // Particle loop

    while ((p = dynamic_cast<GHepParticle *>(iter.Next())))
    {
      int myStatus = p->Status();
      int myPdg = p->Pdg();
      double myCharge = p->Charge();
      double myE = p->E();
      double myKE = p->KinE();

      if (myStatus == 1 && myCharge == 0)
      {
        if (myPdg == kPdgNeutron)
        {
          e0Sum += myKE;
        }
        else
        {
          e0Sum += myE;
        }
      }
    }

    // Fill the bins

    myHist->Fill(e0Sum);
    mcrec->Clear();
  }

  // Normalization and stat box removal

  TH1 *normh = (TH1D *)(myHist->Clone("normh"));
  normh->Scale(1. / normh->GetEntries());
  normh->SetStats(false);

  // Create the output file

  TFile *outfile = new TFile("missing_energy.root", "RECREATE");
  outfile->cd();
  normh->Write();
  outfile->Close();
}
