//____________________________________________________________________________
/*
\brief    Plots a 1D histogram depicting the distribution of the total final
          state transverse momentum of each event, restricted to events with
          at least 1 π+ in the final state.
          
          To see the documentation for this example, visit: 
          https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=349
          
\author   Ishaan Vohra <ivohra@exeter.edu / ishaanklv@gmail.com>
          Phillips Exeter Academy
\created  August 16, 2022
*/
//____________________________________________________________________________

#include <iostream>
#include <string>
#include <iomanip>
#include <map>
#include <random>
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
#include "TVector3.h"

using namespace genie;

void transverse_p_pip_only(string filename = "/hepstore/ivohra/miniboone1.ghep.root")
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

  auto myHist = new TH1D("h1", "Total Transverse Momentum for Events Producing At Least 1 Final State π+;Momentum (GeV/c);Probability", 100, 0, 1);

  // Event loop

  for (Long64_t i = 0; i < (tree->GetEntries()); i++)
  {
    tree->GetEntry(i);

    EventRecord &event = *(mcrec->event);
    TObjArrayIter iter(&event);
    GHepParticle *p = 0;

    double Ptotx = 0;
    double Ptoty = 0;

    // Count the number of π+

    int pipCounter = 0;

    // Particle loop

    while ((p = dynamic_cast<GHepParticle *>(iter.Next())))
    {
      int myStatus = p->Status();
      int myPdg = p->Pdg();
      double myPx = p->Px();
      double myPy = p->Py();

      if (myStatus == 1)
      {
        Ptotx += myPx;
        Ptoty += myPy;

        if (myPdg == 211)
        {
          pipCounter++;
        }
      }
    }

    // Fill the bins

    if (pipCounter)
    {
      TVector3 Ptotvec(Ptotx, Ptoty, 0);
      myHist->Fill(Ptotvec.Mag());
      mcrec->Clear();
    }
  }

  // Normalization and stat box removal

  TH1 *normh = (TH1D *)(myHist->Clone("normh"));
  normh->Scale(1. / normh->GetEntries());
  normh->SetStats(false);

  // Create the output file

  TFile *outfile = new TFile("transverse_p_pip_only.root", "RECREATE");
  outfile->cd();
  normh->Write();
  outfile->Close();
}
