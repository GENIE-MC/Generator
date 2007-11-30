//____________________________________________________________________________
/*!

\class    genie::GViewerMainFrame

\brief    GENIE Viewer Main Frame

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 07, 2004

*/
//____________________________________________________________________________

#ifndef _GVIEWER_MAIN_FRAME_H_
#define _GVIEWER_MAIN_FRAME_H_

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
class TF1;

namespace genie {

class GEVGDriver;
class EventRecord;
class GHepPrinter;
class GHepDrawer;

class GViewerMainFrame : public TGMainFrame {

public:

   GViewerMainFrame(const TGWindow * p, UInt_t w, UInt_t h);
   virtual ~GViewerMainFrame();

   void Close              (void) { gApplication->Terminate(0); }
   void Exit               (void) { Close();                    }
   void NextEvent          (void);
   void ShowEvent          (EventRecord * ev_rec);

private:

   void                Initialize                (void);
   void                DefineLayoutHints         (void);
   TGTab *             BuildViewerTabs           (void);
   TGGroupFrame *      BuildGenieControls        (void);
   TGGroupFrame *      BuildImageButtonFrame     (void);
   TGGroupFrame *      BuildNeutrinoP4Controls   (void);
   TGGroupFrame *      BuildInitialStateControls (void);
   const char *        Icon         (const char * name);

   //-- GEVGDriver driver

   GEVGDriver * fEVGDriver;

   //-- GUI widgets

   TGMainFrame *            _main;
   TGGroupFrame *           fControlsGroupFrame;
   TGGroupFrame *           fNuP4ControlsGroupFrame;
   TGGroupFrame *           fInitStateGroupFrame;
   TGGroupFrame *           fImgButtonGroupFrame;
   TGCompositeFrame *       fMainFrame;
   TGCompositeFrame *       fUpperFrame;
   TGCompositeFrame *       fLowerFrame;
   TGCompositeFrame *       fLowerLeftFrame;
   TGCompositeFrame *       fLowerRightFrame;
   TGCompositeFrame *       fFeynmanTab;
   TGCompositeFrame *       fGHepTab;
   TGTab *                  fViewerTabs;
   TRootEmbeddedCanvas *    fEmbeddedCanvas;
   TGTextEdit *             fGHep;
   TGCheckButton *          fCcCheckButton;
   TGCheckButton *          fNcCheckButton;
   TGCheckButton *          fQelCheckButton;
   TGCheckButton *          fSppCheckButton;
   TGCheckButton *          fDisCheckButton;
   TGCheckButton *          fCohCheckButton;
   TGStatusBar *            fStatusBar;
   TGLayoutHints *          fFeynmanTabLayout;
   TGLayoutHints *          fGHepTabLayout;
   TGLayoutHints *          fStatusBarLayout;
   TGLayoutHints *          fLowerLeftFrameLayout;
   TGLayoutHints *          fLowerRightFrameLayout;
   TGLayoutHints *          fControlsFrameLayout;
   TGLayoutHints *          fViewerTabsLayout;
   TGMatrixLayout *         fButtonMatrixLayout;
   TGMatrixLayout *         fNuP4MatrixLayout;
   TGMatrixLayout *         fInitStateMatrixLayout;
   TGPictureButton *        fNextEventButton;
   TGPictureButton *        fExitButton;
   TGNumberEntry *          fPx;
   TGNumberEntry *          fPy;
   TGNumberEntry *          fPz;
   TGNumberEntry *          fE;
   TGNumberEntry *          fA;
   TGNumberEntry *          fZ;
   TGLabel *                fPxLabel;
   TGLabel *                fPyLabel;
   TGLabel *                fPzLabel;
   TGLabel *                fELabel;
   TGLabel *                fALabel;
   TGLabel *                fZLabel;
   TGLabel *                fNuLabel;
   TGLabel *                fEmptyLabel;
   TGLabel *                fEmptyLabel2;
   TGComboBox *             fNu;

   GHepDrawer *             fGHepDrawer;
   GHepPrinter *            fGHepPrinter;

   ClassDef(GViewerMainFrame, 0)
};

}       // genie namespace
#endif  // _GVIEWER_MAIN_FRAME_H_

