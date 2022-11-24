//____________________________________________________________________________
/*
\brief    Plot a normalized 1D histogram depicting the distribution of Pi+ 
          produced in each event
\author   Ishaan Vohra <ishaanklv@gmail.com>
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

#include "TString.h" 
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
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

void pip_distribution(string filename = "/hepstore/ivohra/miniboone1.ghep.root")
{

// Open the GHEP ROOT file

  TFile infile(filename.c_str());


  // Get the tree header


  NtpMCTreeHeader * header =
    dynamic_cast<NtpMCTreeHeader*> (infile.Get("header"));

  // Get the GENIE GHEP tree and set its branch address

  TTree * tree = dynamic_cast<TTree*> (infile.Get("gtree"));
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  // Make histogram

auto myHist = new TH1D("h1","Distribution of Number of π+ Produced in Each Event",6,-0.5,5.5);

  // Event loop

  for(Long64_t i=0; i < (tree->GetEntries()); i++){
    tree->GetEntry(i);

 EventRecord & event = *(mcrec->event);

   TObjArrayIter iter(&event);
   GHepParticle * p = 0;

    int counter = 0;

   // Particle loop

   while((p = dynamic_cast<GHepParticle *>(iter.Next()))) { 

      int pdgc   = p->Pdg();
      int status = p->Status();

      if(status == 1 && pdgc == 211){

      counter++;

        } 
    }

    myHist->Fill((counter));
    mcrec->Clear();
  }

//normalization

 TH1*normh = (TH1D*)(myHist->Clone("normh"));
   normh->Scale(1./normh->GetEntries());
   normh->SetStats(false);


TFile *outfile = new TFile("pip_distribution.root","RECREATE"); 
outfile->cd();
normh->Write();   
outfile->Close();

}


