//____________________________________________________________________________
/*!

\class    genie::GViewerMainFrame

\brief    GENIE Viewer Main Frame

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 07, 2004

\cpright  Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GVIEWER_MAIN_FRAME_H_
#define _GVIEWER_MAIN_FRAME_H_

#include <string>

#include <TApplication.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <RQ_OBJECT.h>

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

using std::string;

namespace genie {

class EventRecord;
class NtpMCEventRecord;

namespace gview {

class MCTruthDisplay;

class GViewerMainFrame : public TGMainFrame {

public:
   GViewerMainFrame(const TGWindow * p, UInt_t w, UInt_t h);
   virtual ~GViewerMainFrame();

   void Close     (void) { gApplication->Terminate(0); }
   void Exit      (void) { Close();                    }
   void FileOpen  (void);
   void NextEvent (void);
   void ShowEvent (EventRecord * ev_rec);

private:

   void Init(void);
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

   ClassDef(GViewerMainFrame, 0)
};

}  // gview namespace
}  // genie namespace

#endif  // _GVIEWER_MAIN_FRAME_H_

