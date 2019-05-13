//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TVirtualX.h>
#include <TGListBox.h>
#include <TGComboBox.h>
#include <TGClient.h>
#include <TGIcon.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGTextEntry.h>
#include <TGMsgBox.h>
#include <TGMenu.h>
#include <TGCanvas.h>
#include <TGTab.h>
#include <TGFileDialog.h>
#include <TGTextEdit.h>
#include <TGStatusBar.h>
#include <TGProgressBar.h>
#include <TGColorSelect.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TRootEmbeddedCanvas.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TLatex.h>
#include <TStyle.h>

#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Tools/Masterclass/GNuMcMainFrame.h"
#include "Tools/Masterclass/MCTruthDisplay.h"

using std::ostringstream;
using std::setprecision;
using std::string;
using std::vector;

using namespace genie;
using namespace genie::masterclass;

ClassImp(GNuMcMainFrame)

//______________________________________________________________________________
GNuMcMainFrame::GNuMcMainFrame(const TGWindow * p, UInt_t w, UInt_t h) :
TGMainFrame(p, w, h)
{
  this->Init();
  this->BuildGUI(p,w,h);
  this->BuildHelpers();
}
//______________________________________________________________________________
void GNuMcMainFrame::Init(void)
{
   fMain                = 0;
   fImgButtonGroupFrame = 0;
   fMainFrame           = 0;
   fUpperFrame          = 0;
   fLowerFrame          = 0;
   fViewerTabs          = 0;
   fFeynmanTab          = 0;
   fGHepTab             = 0;
   fEmbeddedCanvas      = 0;
   fGHep                = 0;
   fStatusBar           = 0;
   fFeynmanTabLayout    = 0;
   fGHepTabLayout       = 0;
   fStatusBarLayout     = 0;
   fViewerTabsLayout    = 0;
   fButtonMatrixLayout  = 0;
   fFileOpenButton      = 0;
   fNextEventButton     = 0;
   fExitButton          = 0;
   fViewTabWidth        = 0;
   fViewTabHeight       = 0;
   
   fTruthDisplay = 0;
   
   fEventFilename = "";
   fEventFile     = 0;
   fGHepTree      = 0;
   fMCRecord      = 0;
   fNuOfEvents    = 0;
   fCurrEventNu   = 0;

}
//______________________________________________________________________________
void GNuMcMainFrame::BuildGUI(const TGWindow * p, UInt_t w, UInt_t h)
{
  fMain = new TGMainFrame(p,w,h);

  fMain->Connect(
     "CloseWindow()", "genie::GNuMcMainFrame", this, "Close()");

  this->BuildMainFrames();

  //
  // UPPER FRAME: add image buttons frame
  //

  fImgButtonGroupFrame = this->BuildImageButtonFrame();
  fUpperFrame -> AddFrame( fImgButtonGroupFrame );

  this->BuildTabs();
  this->BuildStatusBar();

  // initialize
  fMain->SetWindowName("GENIE Event Viewer");
  fMain->MapSubwindows();
  fMain->Resize( fMain->GetDefaultSize() );
  fMain->MapWindow();
}
//______________________________________________________________________________
GNuMcMainFrame::~GNuMcMainFrame()
{
  fMain->Cleanup();
  delete fMain;

  delete fTruthDisplay;
}
//______________________________________________________________________________
void GNuMcMainFrame::BuildMainFrames(void)
{
  fMainFrame  = new TGCompositeFrame(fMain,      1, 1, kVerticalFrame  );
  fUpperFrame = new TGCompositeFrame(fMainFrame, 3, 3, kHorizontalFrame);
  fLowerFrame = new TGCompositeFrame(fMainFrame, 3, 3, kHorizontalFrame);

  fMainFrame -> AddFrame ( fUpperFrame );
  fMainFrame -> AddFrame ( fLowerFrame );
  fMain      -> AddFrame ( fMainFrame  );
}
//______________________________________________________________________________
TGGroupFrame * GNuMcMainFrame::BuildImageButtonFrame(void)
{
  TGGroupFrame * bf = new TGGroupFrame(
        fUpperFrame, "Viewer Control Buttons", kHorizontalFrame);

  fFileOpenButton  = 
    new TGPictureButton(bf, gClient->GetPicture(Icon("open"),32,32));
  fNextEventButton = 
    new TGPictureButton(bf, gClient->GetPicture(Icon("next"),32,32));
  fExitButton = 
    new TGPictureButton(bf, gClient->GetPicture(Icon("exit"), 32,32),
    "gApplication->Terminate(0)");

  fFileOpenButton   -> SetToolTipText( "Open event file" , 1);
  fNextEventButton  -> SetToolTipText( "Get next event" ,  1);
  fExitButton       -> SetToolTipText( "Exit",             1);

  fFileOpenButton  -> Connect(
     "Clicked()","genie::masterclass::GNuMcMainFrame", this,"FileOpen()");
  fNextEventButton  -> Connect(
     "Clicked()","genie::masterclass::GNuMcMainFrame", this,"NextEvent()");

  bf -> AddFrame( fFileOpenButton   );
  bf -> AddFrame( fNextEventButton  );
  bf -> AddFrame( fExitButton       );

  return bf;
}
//______________________________________________________________________________
void GNuMcMainFrame::BuildTabs(void)
{
  fViewTabWidth  = 780;
  fViewTabHeight = 300;

  fViewerTabs = new TGTab(fLowerFrame, 1, 1);

  this->BuildMCTruthTab          ();
  this->BuildFastSimScintCaloTab ();
  this->BuildFastSimCherenkovTab ();

  ULong_t hintViewerTabsLayout = 
       kLHintsTop | kLHintsExpandX | kLHintsExpandY;
  fViewerTabsLayout      
       = new TGLayoutHints(hintViewerTabsLayout, 5, 5, 10, 1);

  fLowerFrame -> AddFrame ( fViewerTabs, fViewerTabsLayout );
}
//______________________________________________________________________________
void GNuMcMainFrame::BuildMCTruthTab(void)
{
// Add tab for displaying MC truth
//
  TGCompositeFrame * tf = 0;

  unsigned int w = fViewTabWidth;
  unsigned int h = fViewTabHeight;

  // tab: Draw "Feynman" diagram

  tf = fViewerTabs->AddTab( "Feynman Diagram" );

  fFeynmanTab = new TGCompositeFrame(tf, w, h, kVerticalFrame);
  fEmbeddedCanvas =  new TRootEmbeddedCanvas("fEmbeddedCanvas", fFeynmanTab, w, h);

  fEmbeddedCanvas -> GetCanvas() -> SetBorderMode (0);
  fEmbeddedCanvas -> GetCanvas() -> SetFillColor  (0);

  ULong_t hintFeynmanTabLayout  = 
     kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY;
  fFeynmanTabLayout = 
     new TGLayoutHints(hintFeynmanTabLayout, 5, 5, 10, 1);

  fFeynmanTab -> AddFrame( fEmbeddedCanvas, fFeynmanTabLayout );
  tf          -> AddFrame( fFeynmanTab,     fFeynmanTabLayout );

  // tab: Print GHEP record

  tf = fViewerTabs->AddTab("GHEP Record");

  fGHepTab = new TGCompositeFrame(tf, w, h, kVerticalFrame);

  fGHep = new TGTextEdit(fGHepTab, w, h, kSunkenFrame | kDoubleBorder);
  fGHep->AddLine( "GHEP:" );

  ULong_t hintGHepTabLayout = 
     kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY;
  fGHepTabLayout         
     = new TGLayoutHints(hintGHepTabLayout, 5, 5, 10, 1);

  fGHepTab -> AddFrame(fGHep,    fGHepTabLayout);
  tf       -> AddFrame(fGHepTab, fGHepTabLayout);
}
//______________________________________________________________________________
void GNuMcMainFrame::BuildFastSimScintCaloTab (void)
{
// Build tab for displaying fast simulation results of scintillator calorimeter
// response and for controlling simulation inputs.
//
  TGCompositeFrame * tf = 0;

  tf = fViewerTabs->AddTab("FastSim/ScintCalo");

}
//______________________________________________________________________________
void GNuMcMainFrame::BuildFastSimCherenkovTab (void)
{
// Build tab for displaying fast simulation results of Cherenkov detector 
// response and for controlling simulation inputs.
//
  TGCompositeFrame * tf = 0;

  tf = fViewerTabs->AddTab("FastSim/Cherenkov");

}
//______________________________________________________________________________
void GNuMcMainFrame::BuildStatusBar(void)
{
  Int_t parts[] = { 60, 20, 20 };
  fStatusBar = new TGStatusBar(fMain, 50, 10, kHorizontalFrame);
  fStatusBar->SetParts(parts, 3);

  ULong_t hintStatusBarLayout = 
     kLHintsBottom | kLHintsLeft | kLHintsExpandX;
  fStatusBarLayout       
     = new TGLayoutHints(hintStatusBarLayout, 0, 0,  2, 0);

  fMain->AddFrame(fStatusBar, fStatusBarLayout);
}
//______________________________________________________________________________
const char * GNuMcMainFrame::Icon(const char * name)
{
  ostringstream pic;
  pic  << gSystem->Getenv("GENIE") << "/data/icons/" << name << ".xpm";

  LOG("MasterClass", pINFO) << "Loading icon: " << pic.str();

  return pic.str().c_str();
}
//______________________________________________________________________________
void GNuMcMainFrame::BuildHelpers(void)
{
  fTruthDisplay = new MCTruthDisplay(fEmbeddedCanvas,fGHep);
}
//______________________________________________________________________________
void GNuMcMainFrame::FileOpen(void)
{
  fStatusBar->SetText( "Asking for event file name...", 0);

  static TString dir(".");
  const char * kFileExt[] = {"GHEP/ROOT event files", "*.root", 0, 0};

  TGFileInfo fi;
  fi.fFileTypes = kFileExt;
  fi.fIniDir    = StrDup(dir.Data());

  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);

  if( fi.fFilename ) {
     fEventFilename = string( fi.fFilename );

     ostringstream cmd;
     cmd << "Will read events from: " << fEventFilename;
     fStatusBar -> SetText( cmd.str().c_str(), 0 );

     if(fEventFile) {
       fEventFile->Close();
       delete fEventFile;
     }
     if(fGHepTree) {
       delete fGHepTree;
     }

     fEventFile = 
         new TFile(fEventFilename.c_str(),"READ");
     fGHepTree = 
         dynamic_cast <TTree *> (fEventFile->Get("gtree"));
     if(!fGHepTree) {
        LOG("MasterClass", pFATAL) 
            << "No GHEP event tree in input file: " << fEventFilename;
        gAbortingInErr=true;
        exit(1);
     }
     fCurrEventNu = 0;
     fNuOfEvents  = fGHepTree->GetEntries();
     LOG("MasterClass", pNOTICE)  
       << "Input GHEP event tree has " << fNuOfEvents 
       << ((fNuOfEvents==1) ? " entry." : " entries.");

     NtpMCTreeHeader * thdr = 
         dynamic_cast <NtpMCTreeHeader *> ( fEventFile->Get("header") );
     LOG("MasterClass", pNOTICE) 
         << "Input tree header: " << *thdr;

     fGHepTree->SetBranchAddress("gmcrec", &fMCRecord);

  }
}
//______________________________________________________________________________
void GNuMcMainFrame::NextEvent(void)
{
  if(fCurrEventNu >= fNuOfEvents-1) {
	exit(1);
  }

  fGHepTree->GetEntry(fCurrEventNu);
  fCurrEventNu++;

  EventRecord * event = fMCRecord->event;

  this->ShowEvent(event);
}
//______________________________________________________________________________
void GNuMcMainFrame::ShowEvent(EventRecord * event)
{
  fTruthDisplay->DrawDiagram(event);
  fTruthDisplay->PrintEventRecord(event);
}
//______________________________________________________________________________
