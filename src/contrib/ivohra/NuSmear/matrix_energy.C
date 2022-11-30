//____________________________________________________________________________
/*
\brief    Plots a 2D histogram depicting the energy smearing matrix given by
          the chosen NuSmear preset.
          
          To see the full NuSmear paper, visit:
          https://inspirehep.net/literature/2150455
          
          To see the documentation for this code, visit: 
          https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=349
          
\author   Ishaan Vohra <ivohra@exeter.edu / ishaanklv@gmail.com>
          Phillips Exeter Academy
\created  August 16, 2022
*/
//____________________________________________________________________________

#include <iostream>
#include <string>
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
#include "TStyle.h"
#include "NuSmear.h"

using namespace genie;

void matrix_energy(string filename = "/hepstore/ivohra/miniboone1.ghep.root")
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

  auto myHist = new TH2D("h1", "Energy Smearing Matrix;True Energy (GeV);Reconstructed Energy (GeV)", 100, 0, 7, 100, 0, 7);

  // Event loop

  for (Long64_t i = 0; i < (tree->GetEntries()); i++)
  {
    tree->GetEntry(i);

    EventRecord &event = *(mcrec->event);
    TObjArrayIter iter(&event);
    GHepParticle *p = 0;

    double eSumPreSm = 0;
    double eSumPosSm = 0;

    // Particle loop

    while ((p = dynamic_cast<GHepParticle *>(iter.Next())))
    {
      int myStatus = p->Status();
      int myPdg = p->Pdg();
      double myE = p->E();
      double myPx = p->Px();
      double myPy = p->Py();
      double myPz = p->Pz();
      double myKE = p->KinE();

      if (myStatus == 1 && !(pdg::IsNeutrino(myPdg)) && !(pdg::IsIon(myPdg)))
      {
        if (pdg::IsNeutronOrProton(myPdg))
        {
          eSumPreSm += myKE;
        }
        else
        {
          eSumPreSm += myE;
        }
        eSumPosSm += smearE(myPdg, myE, myKE, myPx, myPy, myPz, "default");
      }
    }

    // Fill the bins

    myHist->Fill(eSumPreSm, eSumPosSm);
    mcrec->Clear();
  }

  // Normalization and stat box removal

  TH2 *normh = (TH2D *)(myHist->Clone("normh"));
  normh->Scale(1. / normh->GetEntries());
  normh->SetStats(false);

  // Create the output file

  TFile *outfile = new TFile("matrix_energy.root", "RECREATE");
  outfile->cd();
  normh->Write();
  outfile->Close();
}
