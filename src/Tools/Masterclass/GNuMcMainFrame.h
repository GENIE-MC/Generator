//____________________________________________________________________________
/*!

\class    genie::GNuMcMainFrame

\brief    GENIE Neutrino Masterclass app main frame

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 07, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_NUMC_MAIN_FRAME_H_
#define _G_NUMC_MAIN_FRAME_H_

#include <string>

#include <TApplication.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <RQ_OBJECT.h>
/*
class TApplication;
class TVirtualX;
class TSystem;
class TGListBox;
class TGComboBox;
class TGClient;
class TGIcon;
class TGLabel;
class TGNumberEntry;
class TGTextEntry;
class TGMsgBox;
class TGMenu;
class TGCanvas;
class TGTab;
class TGFileDialog;
class TGTextEdit;
class TGStatusBar;
class TGProgressBar;
class TGColorSelect;
class TCanvas;
class TGraphAsymmErrors;
class TRootEmbeddedCanvas;
class TFile;
class TTree;
*/

//class TApplication;
#include <TVirtualX.h>
#include <TSystem.h>
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
#include <TFile.h>
#include <TTree.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Tools/Masterclass/MCTruthDisplay.h"

using std::string;

namespace genie {
namespace masterclass {

class GNuMcMainFrame : public TGMainFrame {

public:
   GNuMcMainFrame(const TGWindow * p, UInt_t w, UInt_t h);
   virtual ~GNuMcMainFrame();

   void Close     (void) { gApplication->Terminate(0); }
   void Exit      (void) { Close();                    }
   void FileOpen  (void);
   void NextEvent (void);
   void ShowEvent (EventRecord * ev_rec);

private:

   void           Init                     (void);
   void           BuildHelpers             (void);
   void           BuildGUI                 (const TGWindow * p, UInt_t w, UInt_t h);
   void           BuildMainFrames          (void);
   void           BuildTabs                (void);
   void           BuildMCTruthTab          (void);
   void           BuildFastSimScintCaloTab (void);
   void           BuildFastSimCherenkovTab (void);
   void           BuildStatusBar           (void);
   TGGroupFrame * BuildImageButtonFrame    (void);
   const char *   Icon                     (const char * name);

   // GUI widgets & properties
   TGMainFrame *            fMain;
   TGGroupFrame *           fImgButtonGroupFrame;
   TGCompositeFrame *       fMainFrame;
   TGCompositeFrame *       fUpperFrame;
   TGCompositeFrame *       fLowerFrame;
   TGTab *                  fViewerTabs;
   TGCompositeFrame *       fFeynmanTab;
   TGCompositeFrame *       fGHepTab;
   TRootEmbeddedCanvas *    fEmbeddedCanvas;
   TGTextEdit *             fGHep;
   TGStatusBar *            fStatusBar;
   TGLayoutHints *          fFeynmanTabLayout;
   TGLayoutHints *          fGHepTabLayout;
   TGLayoutHints *          fStatusBarLayout;
   TGLayoutHints *          fViewerTabsLayout;
   TGMatrixLayout *         fButtonMatrixLayout;
   TGPictureButton *        fFileOpenButton;
   TGPictureButton *        fNextEventButton;
   TGPictureButton *        fExitButton;
   unsigned int             fViewTabWidth;
   unsigned int             fViewTabHeight;

   // utility classes
   MCTruthDisplay * fTruthDisplay;

   // input events
   string             fEventFilename;
   TFile*             fEventFile;
   TTree*             fGHepTree;
   NtpMCEventRecord * fMCRecord;
   Long64_t           fNuOfEvents;
   Long64_t           fCurrEventNu;

   ClassDef(GNuMcMainFrame, 1)
};

}  // masterclass namespace
}  // genie namespace

#endif  // _G_NUMC_MAIN_FRAME_H_

