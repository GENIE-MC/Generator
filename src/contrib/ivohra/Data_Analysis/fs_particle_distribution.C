//____________________________________________________________________________
/*
\brief    Plots a normalized 2D histogram depicting the distribution of the
          number of each type of final state particle produced per event.
          
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
#include "GHepParticle.h"
#include "TObjArray.h"
#include "PDGLibrary.h"
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
#include "TDatabasePDG.h"

using namespace genie;

void fs_particle_distribution(string filename = "/hepstore/ivohra/miniboone1.ghep.root")
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

  auto myHist = new TH2D("h1", "Distribution of the Number of Each Type of Final State Particle Produced per Event;Number Produced per Event", 30, -0.5, 29.5, 30, 0.5, 30.5);

  // map PDG to bins

  std::map<int, int> pdgToBin;

  pdgToBin[kPdgPiP] = 1;
  pdgToBin[kPdgPiM] = 2;
  pdgToBin[kPdgPi0] = 3;
  pdgToBin[kPdgEta] = 4;
  pdgToBin[kPdgEtaPrm] = 5;
  pdgToBin[kPdgEtac] = 6;
  pdgToBin[kPdgEtab] = 7;
  pdgToBin[kPdgRhoP] = 8;
  pdgToBin[kPdgRhoM] = 9;
  pdgToBin[kPdgRho0] = 10;
  pdgToBin[kPdgPhi] = 11;
  pdgToBin[kPdgJpsi] = 12;
  pdgToBin[kPdgKP] = 13;
  pdgToBin[kPdgKM] = 14;
  pdgToBin[kPdgK0] = 15;
  pdgToBin[kPdgAntiK0] = 16;
  pdgToBin[kPdgKStarP] = 17;
  pdgToBin[kPdgKStarM] = 18;
  pdgToBin[kPdgKStar0] = 19;
  pdgToBin[kPdgDP] = 20;
  pdgToBin[kPdgDM] = 21;
  pdgToBin[kPdgD0] = 22;
  pdgToBin[kPdgAntiD0] = 23;
  pdgToBin[kPdgDPs] = 24;
  pdgToBin[kPdgDMs] = 25;
  pdgToBin[kPdgGamma] = 26;
  pdgToBin[kPdgProton] = 27;
  pdgToBin[kPdgAntiProton] = 28;
  pdgToBin[kPdgMuon] = 29;
  pdgToBin[kPdgAntiMuon] = 30;

  // Y-axis labelling

  TDatabasePDG *db = TDatabasePDG::Instance();

  for (auto &x : pdgToBin)
  {

    TParticlePDG *part = db->GetParticle(x.first);
    myHist->GetYaxis()->SetBinLabel(x.second, part->GetName());
  }

  // Event loop

  for (Long64_t i = 0; i < (tree->GetEntries()); i++)
  {
    tree->GetEntry(i);

    EventRecord &event = *(mcrec->event);
    TObjArrayIter iter(&event);
    GHepParticle *p = 0;

    int counter[31] = {0};

    // Particle loop

    while ((p = dynamic_cast<GHepParticle *>(iter.Next())))
    {
      int pdgc = p->Pdg();
      int status = p->Status();

      if (status == 1)
      {
        auto binNumber = (pdgToBin.find(pdgc))->second;
        counter[binNumber]++;
      }
    }

    // Fill the bins

    for (int i = 1; i <= 30; i++)
    {
      myHist->Fill(counter[i], i);
    }
    mcrec->Clear();
  }

  // Normalization and stat box removal

  TH2 *normh = (TH2D *)(myHist->Clone("normh"));
  normh->Scale(1. / normh->GetEntries());
  normh->SetStats(false);

  // Create the output file

  TFile *outfile = new TFile("fs_particle_distribution.root", "RECREATE");
  outfile->cd();
  normh->Write();
  outfile->Close();
}
