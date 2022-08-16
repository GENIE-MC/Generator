//____________________________________________________________________________
/*
\brief    Plot a 1D histogram depicting the distribution of the total final 
          state energy of each event
\author   Ishaan Vohra <ishaanklv@gmail.com>
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
#include <iomanip>
#include <map>
#include <random>


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


void total_energy()
{

using namespace genie;

// Open the GHEP/ROOT file

  string filename = "/hepstore/ivohra/miniboone1.ghep.root";
  TFile infile(filename.c_str());


  // Get the tree header & print it


  NtpMCTreeHeader * header =
    dynamic_cast<NtpMCTreeHeader*> (infile.Get("header"));

  // Get the GENIE GHEP tree and set its branch address

  TTree * tree = dynamic_cast<TTree*> (infile.Get("gtree"));
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //Make histogram

auto myHist = new TH1D("h1","Total Final State Energy;Energy (GeV);Probability",100,0,6);

  // Event loop

  for(Long64_t i=0; i < (tree->GetEntries()); i++){
    tree->GetEntry(i);

 EventRecord & event = *(mcrec->event);

   TObjArrayIter iter(&event);
   GHepParticle * p = 0;

    double eSum = 0;

    // Particle loop

   while((p = dynamic_cast<GHepParticle *>(iter.Next()))) {

    int myStatus = p->Status();
    double myE = p->E();
    int myPdg = p->Pdg();
    double myKE = p->KinE();

          if(myStatus == 1){

            if(pdg::IsNeutronOrProton(myPdg) || pdg::IsIon(myPdg)){
              eSum += myKE;
            } else {
              eSum += myE;
            }    
          }
      }   

    myHist->Fill(eSum);

    mcrec->Clear();
 
  }

 TH1*normh = (TH1D*)(myHist->Clone("normh"));
   normh->Scale(1./normh->GetEntries());

normh->SetStats(false);

TFile *outfile = new TFile("total_energy.root","RECREATE"); 
outfile->cd();
normh->Write();   
outfile->Close();

}
