//____________________________________________________________________________
/*!

\class    genie::GenieViewer

\brief    GENIE Viewer Main Frame

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 07, 2004

*/
//____________________________________________________________________________

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

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
#include <TF1.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TLatex.h>
#include <TStyle.h>

#include "EVGDrivers/GENIE.h"
#include "EVGCore/EventRecord.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Viewer/RendererQEL.h"
#include "Viewer/GHepPrinter.h"
#include "Viewer/GenieViewer.h"

using std::ostringstream;
using std::setprecision;
using std::string;
using std::vector;

using namespace genie;

ClassImp(GenieViewer)

const char * neutrinos[] = {
             "nu_e", "nu_mu", "nu_tau", "nu_e_bar", "nu_mu_bar", "nu_tau_bar" };

const int neutrino_pdg[] = { 12, 14, 16, -12, -14, -16 };

//______________________________________________________________________________
GenieViewer::GenieViewer(const TGWindow * p, UInt_t w, UInt_t h) :
TGMainFrame(p, w, h)
{
  _main = new TGMainFrame(p,w,h);

  _main->Connect("CloseWindow()", "genie::GenieViewer", this, "Close()");

  //-- define TGLayoutHints

  DefineLayoutHints();

  //-- instantiate main frames (below menu & above the status bar)

  fMainFrame  = new TGCompositeFrame(_main,      1, 1, kVerticalFrame);
  fUpperFrame = new TGCompositeFrame(fMainFrame, 3, 3, kHorizontalFrame);
  fLowerFrame = new TGCompositeFrame(fMainFrame, 3, 3, kHorizontalFrame);

  //-- UPPER FRAME: add image buttons frame

  fImgButtonGroupFrame = BuildImageButtonFrame();

  fUpperFrame -> AddFrame( fImgButtonGroupFrame );

  //-- LOWER FRAME:
  //         left  : Controls
  //         right : Feynman diagram & STDHEP record viewers

  fLowerLeftFrame   = new TGCompositeFrame(fLowerFrame, 3, 3, kVerticalFrame);
  fLowerRightFrame  = new TGCompositeFrame(fLowerFrame, 3, 3, kVerticalFrame);

  fLowerFrame  -> AddFrame ( fLowerLeftFrame,   fLowerLeftFrameLayout  );
  fLowerFrame  -> AddFrame ( fLowerRightFrame,  fLowerRightFrameLayout );

  fControlsGroupFrame     = BuildGenieControls();
  fInitStateGroupFrame    = BuildInitialStateControls();
  fNuP4ControlsGroupFrame = BuildNeutrinoP4Controls();

  fViewerTabs             = BuildViewerTabs();

  fLowerLeftFrame  -> AddFrame ( fControlsGroupFrame,  fControlsFrameLayout  );
  fLowerLeftFrame  -> AddFrame ( fInitStateGroupFrame    );
  fLowerLeftFrame  -> AddFrame ( fNuP4ControlsGroupFrame );
  fLowerRightFrame -> AddFrame ( fViewerTabs,          fViewerTabsLayout     );

  //-- add top/bottom main frames to main frame, and main frame to main

  fMainFrame        -> AddFrame ( fUpperFrame    );
  fMainFrame        -> AddFrame ( fLowerFrame );

  _main             -> AddFrame ( fMainFrame  );

  //-- add Status Bar

  Int_t parts[] = { 60, 20, 20 };

  fStatusBar = new TGStatusBar(_main, 50, 10, kHorizontalFrame);

  fStatusBar->SetParts(parts, 3);

  _main->AddFrame(fStatusBar, fStatusBarLayout);

  //-- initialize

  Initialize();

  _main->SetWindowName("GENIE Viewer");

  _main->MapSubwindows();

  _main->Resize( _main->GetDefaultSize() );

  _main->MapWindow();
}
//______________________________________________________________________________
GenieViewer::~GenieViewer()
{
  _main->Cleanup();
  
  delete _main;
}
//______________________________________________________________________________
void GenieViewer::DefineLayoutHints(void)
{
  ULong_t hintFeynmanTabLayout  = kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintGHepTabLayout   = kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintStatusBarLayout   = kLHintsBottom | kLHintsLeft | kLHintsExpandX;
  ULong_t hintLLeftFrameLayout  = kLHintsCenterY;
  ULong_t hintLRightFrameLayout = kLHintsTop | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintControlsLayout    = kLHintsTop | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintViewerTabsLayout  = kLHintsTop | kLHintsExpandX | kLHintsExpandY;

  fFeynmanTabLayout      = new TGLayoutHints(hintFeynmanTabLayout,    5, 5, 10, 1);
  fGHepTabLayout       = new TGLayoutHints(hintGHepTabLayout,     5, 5, 10, 1);
  fStatusBarLayout       = new TGLayoutHints(hintStatusBarLayout,     0, 0,  2, 0);
  fLowerLeftFrameLayout  = new TGLayoutHints(hintLLeftFrameLayout,    1, 1,  1, 1);
  fLowerRightFrameLayout = new TGLayoutHints(hintLRightFrameLayout,   1, 1,  1, 1);
  fControlsFrameLayout   = new TGLayoutHints(hintControlsLayout,      5, 5, 10, 1);
  fViewerTabsLayout      = new TGLayoutHints(hintViewerTabsLayout,    5, 5, 10, 1);
}
//______________________________________________________________________________
TGGroupFrame * GenieViewer::BuildImageButtonFrame(void)
{
  TGGroupFrame * button_group_frame = new TGGroupFrame(
                       fUpperFrame, "Viewer Control Buttons", kHorizontalFrame);

//  fButtonMatrixLayout = new TGMatrixLayout(button_group_frame, 0, 2, 1);

//  button_group_frame->SetLayoutManager( fButtonMatrixLayout );

  fNextEventButton = new TGPictureButton(button_group_frame,
                                     gClient->GetPicture(Icon("next"),32,32));
  fExitButton      = new TGPictureButton(button_group_frame,
                                     gClient->GetPicture(Icon("exit"),32,32),
                                                 "gApplication->Terminate(0)");

  fNextEventButton  -> SetToolTipText( "Next Event" ,   1);
  fExitButton       -> SetToolTipText( "Exit Viewer",   1);

  fNextEventButton  -> Connect("Clicked()","genie::GenieViewer", this,"NextEvent()");
  
  button_group_frame -> AddFrame( fNextEventButton  );
  button_group_frame -> AddFrame( fExitButton       );

  return button_group_frame;
}
//______________________________________________________________________________
TGGroupFrame * GenieViewer::BuildGenieControls(void)
{
  TGGroupFrame *  controls_frame = new TGGroupFrame(fLowerLeftFrame,
                                               "Run Switches",  kVerticalFrame);

  fCcCheckButton  = new TGCheckButton(controls_frame, "Toggle CC",  101);
  fNcCheckButton  = new TGCheckButton(controls_frame, "Toggle NC",  102);
  fQelCheckButton = new TGCheckButton(controls_frame, "Toggle QEL", 103);
  fSppCheckButton = new TGCheckButton(controls_frame, "Toggle SPP", 104);
  fDisCheckButton = new TGCheckButton(controls_frame, "Toggle DIS", 105);
  fCohCheckButton = new TGCheckButton(controls_frame, "Toggle COH", 106);

  controls_frame -> AddFrame( fCcCheckButton  );
  controls_frame -> AddFrame( fNcCheckButton  );
  controls_frame -> AddFrame( fQelCheckButton );
  controls_frame -> AddFrame( fSppCheckButton );
  controls_frame -> AddFrame( fDisCheckButton );
  controls_frame -> AddFrame( fCohCheckButton );

  return controls_frame;
}
//______________________________________________________________________________
TGGroupFrame * GenieViewer::BuildNeutrinoP4Controls(void)
{
  TGGroupFrame *  nup4_frame = new TGGroupFrame(fLowerLeftFrame,
                                               "Neutrino 4-P",  kVerticalFrame);

  fNuP4MatrixLayout = new TGMatrixLayout(nup4_frame, 0, 2, 5);

  nup4_frame->SetLayoutManager( fNuP4MatrixLayout );

  fPx = new TGNumberEntry(nup4_frame, 0, 7, 3, TGNumberFormat::kNESReal);
  fPy = new TGNumberEntry(nup4_frame, 0, 7, 3, TGNumberFormat::kNESReal);
  fPz = new TGNumberEntry(nup4_frame, 0, 7, 3, TGNumberFormat::kNESReal);
  fE  = new TGNumberEntry(nup4_frame, 0, 7, 3, TGNumberFormat::kNESReal);
                   
  fPxLabel = new TGLabel(nup4_frame, new TGString( "Px = "));
  fPyLabel = new TGLabel(nup4_frame, new TGString( "Py = "));
  fPzLabel = new TGLabel(nup4_frame, new TGString( "Pz = "));
  fELabel  = new TGLabel(nup4_frame, new TGString( "E  = "));
  
  nup4_frame -> AddFrame ( fPxLabel );
  nup4_frame -> AddFrame ( fPx      );
  nup4_frame -> AddFrame ( fPyLabel );
  nup4_frame -> AddFrame ( fPy      );
  nup4_frame -> AddFrame ( fPzLabel );
  nup4_frame -> AddFrame ( fPz      );
  nup4_frame -> AddFrame ( fELabel  );
  nup4_frame -> AddFrame ( fE       );

  fEmptyLabel  = new TGLabel(nup4_frame, new TGString( " "));

  nup4_frame -> AddFrame ( fEmptyLabel  );
  nup4_frame -> AddFrame ( fEmptyLabel  );

  return nup4_frame;
}
//______________________________________________________________________________
TGGroupFrame * GenieViewer::BuildInitialStateControls(void)
{
  TGGroupFrame *  init_state_frame = new TGGroupFrame(fLowerLeftFrame,
                                             "Initial State",  kVerticalFrame);

  fInitStateMatrixLayout = new TGMatrixLayout(init_state_frame, 0, 2, 4);

  init_state_frame->SetLayoutManager( fInitStateMatrixLayout );

  fA = new TGNumberEntry(init_state_frame, 1, 4, 0, TGNumberFormat::kNESInteger);
  fZ = new TGNumberEntry(init_state_frame, 1, 4, 0, TGNumberFormat::kNESInteger);
  
  fALabel  = new TGLabel(init_state_frame, new TGString( "A = "));
  fZLabel  = new TGLabel(init_state_frame, new TGString( "Z = "));
  fNuLabel = new TGLabel(init_state_frame, new TGString( "Nu "));

  fNu = new TGComboBox(init_state_frame, 201);
  
  int i = 0;
  while( neutrinos[i] ) {
    fNu->AddEntry(neutrinos[i], i);
    i++;
  }

  fNu -> Resize (80, 20);

  init_state_frame -> AddFrame ( fALabel  );
  init_state_frame -> AddFrame ( fA       );
  init_state_frame -> AddFrame ( fZLabel  );
  init_state_frame -> AddFrame ( fZ       );
  init_state_frame -> AddFrame ( fNuLabel );
  init_state_frame -> AddFrame ( fNu      );

  fEmptyLabel2  = new TGLabel(init_state_frame, new TGString(" "));

  init_state_frame -> AddFrame ( fEmptyLabel2 );
  init_state_frame -> AddFrame ( fEmptyLabel2 );

  return init_state_frame;
}
//______________________________________________________________________________
TGTab * GenieViewer::BuildViewerTabs(void)
{
  TGCompositeFrame * tf = 0;

  unsigned int width  = 780;
  unsigned int height = 300;

  TGTab * tab = new TGTab(fLowerRightFrame, 1, 1);

  //--- tab: "Plotter"

  tf = tab->AddTab( "Feynman Diagram" );

  fFeynmanTab     = new TGCompositeFrame    (tf, width, height, kVerticalFrame);
  fEmbeddedCanvas = new TRootEmbeddedCanvas (
                                 "fEmbeddedCanvas", fFeynmanTab, width, height);

  fEmbeddedCanvas -> GetCanvas() -> SetBorderMode (0);
  fEmbeddedCanvas -> GetCanvas() -> SetFillColor  (0);

  fFeynmanTab -> AddFrame( fEmbeddedCanvas, fFeynmanTabLayout );
  tf          -> AddFrame( fFeynmanTab,     fFeynmanTabLayout );

  //--- tab "Data Viewer"

  tf = tab->AddTab("GHEP Record");

  fGHepTab = new TGCompositeFrame(tf,  width, height, kVerticalFrame);

  fGHep = new TGTextEdit(fGHepTab,  width, height, kSunkenFrame | kDoubleBorder);
  fGHep->AddLine( "GHEP:" );

  fGHepTab -> AddFrame( fGHep,    fGHepTabLayout);
  tf       -> AddFrame( fGHepTab, fGHepTabLayout);

  return tab;
}
//______________________________________________________________________________
const char * GenieViewer::Icon(const char * name)
{
  ostringstream pic;
  pic  << gSystem->Getenv("GENIE") << "/data/icons/" << name << ".xpm";

  return pic.str().c_str();
}
//______________________________________________________________________________
void GenieViewer::Initialize(void)
{
  //-- build GENIE interface object
  
  fGENIE = new GENIE();

  //-- build Feynman diagram renderers and STDHEP record printer

  fQELRenderer   = new RendererQEL;
  fGHepPrinter = new GHepPrinter;

  fQELRenderer->SetEmbeddedCanvas(fEmbeddedCanvas);
  fGHepPrinter->SetTextEdit(fGHep);
}
//______________________________________________________________________________
void GenieViewer::NextEvent(void)
{
  fGENIE->SetInitialState(14,(int) fZ->GetNumber(), (int) fA->GetNumber());  // Z,A

  double px = fPx->GetNumber();
  double py = fPy->GetNumber();
  double pz = fPz->GetNumber();
  double E  = fE ->GetNumber();

  TLorentzVector nu_p4(px,py,pz,E); // px,py,pz,E (GeV)

  EventRecord * ev_rec = fGENIE->GenerateEvent(nu_p4);

  LOG("gviewer", pINFO) << *ev_rec;

  this->ShowEvent(ev_rec);
}
//______________________________________________________________________________
void GenieViewer::ShowEvent(EventRecord * ev_rec) 
{
  LOG("gviewer", pINFO) << "Drawing Feynman Diagram";

  fGHepPrinter->Print(ev_rec);

  Interaction * interaction = ev_rec->GetInteraction();

  if (interaction->GetProcessInfo().IsQuasiElastic())
  {
    fQELRenderer->DrawDiagram(ev_rec);
  }
  else if (interaction->GetProcessInfo().IsDeepInelastic())
  {

  }
  else if (interaction->GetProcessInfo().IsResonant())
  {

  }
  else 
  {

  }
}
//______________________________________________________________________________
