//____________________________________________________________________________
/*!

\program gtestEventLoop

\brief   Example event loop. Use that as a template for your analysis code.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created May 4, 2004

\cpright Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>
#include <TH1.h>
#include <TText.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "NucleonDecay/NucleonDecayUtils.h"
#include "NucleonDecay/NucleonDecayMode.h"


using std::string;
using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);

int    gOptNEvt;
string gOptInpFilename;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFilename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  if(!tree) return 1;

  // Output histogram file
  std::string filename = std::string(file.GetName());
  std::string histofilename = filename.substr(0,filename.size()-5) +
    ".histo.root";
  TFile *histofile = new TFile(histofilename.c_str(),"RECREATE");
 
  // Decayed nucleon histograms
  TH1D* dnPdgHisto = new TH1D("DecayedNucleonPdg", "DecayedNucleonPdg", 5000, -2500, 2500);

  TH1D* dnPHisto = new TH1D("DecayedNucleonMomentum", "DecayedNucleonMomentum [GeV/c]", 100, 0., 0.5);  // [GeV/c]

  TH1D* dnRemovalEnergyHisto = new TH1D("DecayedNucleonRemovalEnergy", "DecayedNucleonRemovalEnergy [GeV]", 100, 0., 0.05); // [GeV/c]

  // Histograms for decayed nucleon daughters
  TH1D* dndaughtersNHisto = new TH1D("DecayedNucleonNDaughters", "DecayedNucleonNDaughters", 6, 0, 6);

  TH1D* dndaughter0PdgHisto = new TH1D("DecayedNucleonDaughter0Pdg", "DecayedNucleonDaughter0Pdg", 1000, -500, 500);
  TH1D* dndaughter1PdgHisto = new TH1D("DecayedNucleonDaughter1Pdg", "DecayedNucleonDaughter1Pdg", 1000, -500, 500);
  TH1D* dndaughter2PdgHisto = new TH1D("DecayedNucleonDaughter2Pdg", "DecayedNucleonDaughter2Pdg", 1000, -500, 500);

    TH1D* dndaughter0PHisto = new TH1D("DecayedNucleonDaughter0Momentum", "DecayedNucleonDaughter0Momentum [GeV/c]", 100, 0., 1.);  // [GeV/c]
    TH1D* dndaughter1PHisto = new TH1D("DecayedNucleonDaughter1Momentum", "DecayedNucleonDaughter1Momentum [GeV/c]", 100, 0., 1.);  // [GeV/c]
    TH1D* dndaughter2PHisto = new TH1D("DecayedNucleonDaughter2Momentum", "DecayedNucleonDaughter2Momentum [GeV/c]", 100, 0., 1.);  // [GeV/c]

    // Histograms for stable final state particles
    TH1D* finalparticlesNHisto = new TH1D("NFinalParticles", "NFinalParticles", 20, 0, 20);
  TH1D* finalparticlesPdgHisto = new TH1D("FinalParticlesPdg", "FinalParticlesPdg", 5000, -2500, 2500);
  TH1D* finalparticlesPHisto = new TH1D("FinalParticlesMomentum", "FinalParticlesMomentum [GeV/c]", 100, 0., 1.);  // [GeV/c]    


  // address for input file
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (gOptNEvt > 0) ?
        TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  //
  // Loop over all events
  //
  // boolean to set, only for the first event, a few TTexts to be written into the histogram filename
  bool first = true;
  TText decayName = TText();
  decayName.SetName("DecayName");
  TText targetName = TText();
  targetName.SetName("TargetName");

  // decayed nucleon index. In the event record, the index 0 element is typically the target, while the index 1 element is the decayed nucleon within the target
  int dnIndex = 1;

  int ndaughters;
  int nfinalparticles;

  for(int i = 0; i < nev; i++) {

    // get next tree entry
    tree->GetEntry(i);

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);

    LOG("testNucleonDecay", pNOTICE) << event;

    if (first) {
      first = !first;

      NucleonDecayMode_t ndm =  (NucleonDecayMode_t)event.Summary()->ExclTagPtr()->DecayMode();
      int npdg =  event.Summary()->InitStatePtr()->TgtPtr()->HitNucPdg();
      string decayNameString = genie::utils::nucleon_decay::AsString(ndm,npdg);
      decayName.SetTitle(decayNameString.c_str());
      string targetNameString = event.Summary()->InitStatePtr()->TgtPtr()->AsString(); 
      targetName.SetTitle(targetNameString.c_str());
    }

    GHepParticle * p = 0;
    TIter event_iter(&event);

    // decayed nucleon histograms
    GHepParticle* dnpart = event.Particle(dnIndex);

    if (dnpart->Status() != 3) {
      LOG("testNucleonDecay", pFATAL) << "Unexpected status code for candidate decayed nucleon: " << dnpart->Status() << ", exiting.";
      exit(1);
    }

    if (dnpart->FirstMother() != 0) {
      LOG("testNucleonDecay", pFATAL) << "Unexpected first mother index for candidate decayed nucleon: " << dnpart->FirstMother() << ", exiting.";
      exit(1);
    }

    dnPdgHisto->Fill(dnpart->Pdg());
    dnPHisto->Fill(dnpart->P4()->Vect().Mag());
    dnRemovalEnergyHisto->Fill(dnpart->RemovalEnergy());

    // histograms for decayed nucleon daughters
    ndaughters = 0;

    while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
    {
      if (p->FirstMother() == dnIndex) {
	if (ndaughters == 0) {
	  dndaughter0PdgHisto->Fill(p->Pdg());
	  dndaughter0PHisto->Fill(p->P4()->Vect().Mag());
	} else if (ndaughters == 1) {
	  dndaughter1PdgHisto->Fill(p->Pdg());
	  dndaughter1PHisto->Fill(p->P4()->Vect().Mag());
	} else if (ndaughters == 2) {
	  dndaughter2PdgHisto->Fill(p->Pdg());
	  dndaughter2PHisto->Fill(p->P4()->Vect().Mag());
	}
	ndaughters++;
      }
    }
    dndaughtersNHisto->Fill(ndaughters);

    // histograms for stable final state particles
    nfinalparticles = 0;
    
    event_iter.Reset();
    while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
      {
	if (p->Status() == kIStStableFinalState ) {
	  finalparticlesPdgHisto->Fill(p->Pdg());
	  finalparticlesPHisto->Fill(p->P4()->Vect().Mag());
	  nfinalparticles++;
	}
      }
    finalparticlesNHisto->Fill(nfinalparticles);


    // clear current mc event record
    mcrec->Clear();

  }//end loop over events

  // close input GHEP event file
  file.Close();

  // write TTexts explicitly. No need to do that for histograms, which are written by default
  decayName.Write();
  targetName.Write();
  
  histofile->Write();
  histofile->Close();

  LOG("testNucleonDecay", pNOTICE)  << "Done!";

  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("testNucleonDecay", pINFO) << "Parsing commad line arguments";

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if( parser.OptionExists('f') ) {
    LOG("testNucleonDecay", pINFO) 
       << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("testNucleonDecay", pFATAL) 
        << "Unspecified input filename - Exiting";
    exit(1);
  }

  // number of events to analyse
  if( parser.OptionExists('n') ) {
    LOG("testNucleonDecay", pINFO) 
      << "Reading number of events to analyze";
    gOptNEvt = parser.ArgAsInt('n');
  } else {
    LOG("testNucleonDecay", pINFO)
      << "Unspecified number of events to analyze - Use all";
    gOptNEvt = -1;
  }
}
//_________________________________________________________________________________
