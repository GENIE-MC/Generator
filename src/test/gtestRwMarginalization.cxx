//____________________________________________________________________________
/*!

\program gtestRwMarginalization

\brief   

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created Oct 19, 2009

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>
#include <TArrayD.h>
#include <TH1D.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "ReWeight/GReWeightMarginalize.h"
#include "ReWeight/GSyst.h"

using std::string;
using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);

int    gOptNEvt;
string gOptInpFilename;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  // open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFilename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );
  if(!tree) return 1;

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (gOptNEvt > 0) ?
        TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  //
  // Loop over all events & calculate weights:
  // Modify the pi absorption fraction by +/- 1 sigma and marginalize
  // all other fates (elastic, inelastic, charge exhange, pi production)
  //

  TArrayD piabs_m1sig_weight(nev);
  TArrayD piabs_p1sig_weight(nev);
  
  for(int i = 0; i < nev; i++) {
    tree->GetEntry(i);
    EventRecord & event = *(mcrec->event);

    LOG("test", pNOTICE)  << " ****** current event nu. = " << i;
  //LOG("test", pNOTICE)  << event;

    double piabsm = 
        rew::margin::MarginalizeFates(
             -1,rew::kINukeTwkDial_FrAbs_pi,&event);
    double piabsp = 
        rew::margin::MarginalizeFates(
             +1,rew::kINukeTwkDial_FrAbs_pi,&event);

    LOG("test", pNOTICE)  << " weight (abs + 1 sigma) = " << piabsp;	
    LOG("test", pNOTICE)  << " weight (abs - 1 sigma) = " << piabsm;

    piabs_m1sig_weight.SetAt(piabsm, i);
    piabs_p1sig_weight.SetAt(piabsp, i);

    // clear current mc event record
    mcrec->Clear();

  }//end loop over events


  //
  // Now plot the effect on the muon energy distribution
  // for events with 1 mu and 0 pions with E > some threshold in the final state.
  //

  int    nbins    = 40;
  double Emu_min  =  0.;
  double Emu_max  =  1.;

  TH1D * hEmu_nominal   = new TH1D("hEmu_nominal",  "nominal",                 nbins,Emu_min,Emu_max);
  TH1D * hEmu_abs_m1sig = new TH1D("hEmu_abs_m1sig","pi absorption: -1 sigma", nbins,Emu_min,Emu_max);
  TH1D * hEmu_abs_p1sig = new TH1D("hEmu_abs_p1sig","pi absorption: +1 sigma", nbins,Emu_min,Emu_max);

  hEmu_nominal   -> SetDirectory(0);
  hEmu_abs_m1sig -> SetDirectory(0);
  hEmu_abs_p1sig -> SetDirectory(0);

  // Loop over all events
  for(int i = 0; i < nev; i++) {
    tree->GetEntry(i);
    EventRecord & event = *(mcrec->event);

    double curr_piabs_m1sig_weight = piabs_m1sig_weight[i];
    double curr_piabs_p1sig_weight = piabs_p1sig_weight[i];

    double Epi_thr = 0.100; // GeV
    bool   has_mu  = false;
    int    npi     = 0;
    double Emu     = 0;
    TIter event_iter(&event);
    GHepParticle * p = 0;
    while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
    {
      if(p->Status() != kIStStableFinalState) continue;

      if(p->Pdg() == kPdgMuon && !has_mu) {
         has_mu = true;
         Emu = p->Energy();
      }
      if(pdg::IsPion(p->Pdg()) && p->Energy() > Epi_thr) {
         npi++;
      }
    }//particle loop

    bool pass = has_mu && (npi==0);
    if(pass) {
       hEmu_nominal   -> Fill(Emu, 1.);
       hEmu_abs_m1sig -> Fill(Emu, curr_piabs_m1sig_weight);
       hEmu_abs_p1sig -> Fill(Emu, curr_piabs_p1sig_weight);
    }

    // clear current mc event record
    mcrec->Clear();

  }//event loop

  // close input GHEP event file
  file.Close();

  // save the histograms
  TFile outf("./inuke_systematics_marginalization.root","recreate");
  hEmu_nominal   -> Write();
  hEmu_abs_m1sig -> Write();
  hEmu_abs_p1sig -> Write();
  outf.Close();

  LOG("test", pNOTICE)  << "Done!";

  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("test", pINFO) << "Parsing commad line arguments";

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if ( parser.OptionExists('f') ) {
    LOG("test", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("test", pFATAL) 
        << "Unspecified input filename - Exiting";
    exit(1);
  }

  // number of events to analyse
  if ( parser.OptionExists('n') ) {
    LOG("test", pINFO) << "Reading number of events to analyze";
    gOptNEvt = parser.ArgAsInt('n');
  } else {
    LOG("test", pINFO)
      << "Unspecified number of events to analyze - Use all";
    gOptNEvt = -1;
  }
}
//_________________________________________________________________________________
