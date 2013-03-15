//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
 @ Aug 25, 2009 - CA
   Removed redundant versions of ParserUtils.h and ParserStatus.h in favor of
   the ones in $GENIE/Conventions and $GENIE/Utils. Updated code accordingly.

*/
//____________________________________________________________________________ 

#include <cassert>
#include <string>
#include <sstream>
#include <vector>

#include <TStyle.h>
#include <TSystem.h>
#include <TGListBox.h>
#include <TGComboBox.h>
#include <TGClient.h>
#include <TGIcon.h>
#include <TGLabel.h>
#include <TGButton.h>
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
#include <TGLayout.h>
#include <TCanvas.h>
#include <TGDockableFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TSQLRow.h>
#include <TSQLResult.h>
#include <TSQLServer.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TLatex.h>
#include <TF1.h>
#include <TH2F.h>

#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "Utils/StringUtils.h"
#include "Utils/GUIUtils.h"

#include "ValidationTools/NuVld/DBI.h"
#include "ValidationTools/NuVld/DBTable.h"
#include "ValidationTools/NuVld/DBQueryString.h"
#include "ValidationTools/NuVld/SqlUtils.hh"
#include "ValidationTools/NuVld/NuVldMainFrame.h"
#include "ValidationTools/NuVld/NuVldUserData.h"
#include "ValidationTools/NuVld/GuiTablePrinter.h"
#include "ValidationTools/NuVld/GuiTableRenderer.h"
#include "ValidationTools/NuVld/GuiMsgBox.h"
#include "ValidationTools/NuVld/GuiMultiLineMsgBox.h"
#include "ValidationTools/NuVld/GuiYNQuestionBox.h"
#include "ValidationTools/NuVld/GuiTextEntryDialog.h"
#include "ValidationTools/NuVld/DBConnection.h"
#include "ValidationTools/NuVld/DBConnectionDialog.h"
#include "ValidationTools/NuVld/GuiSysLogSingleton.h"
#include "ValidationTools/NuVld/GuiBrowserSingleton.h"
#include "ValidationTools/NuVld/GraphUtils.hh"
#include "ValidationTools/NuVld/GuiNuDataSelectionTab.h"
#include "ValidationTools/NuVld/GuiElDataSelectionTab.h"
#include "ValidationTools/NuVld/GuiSFDataSelectionTab.h"
#include "ValidationTools/NuVld/NuVldConstants.h"
#include "ValidationTools/NuVld/GuiHelpHandler.h"
#include "ValidationTools/NuVld/GuiStackHandler.h"
#include "ValidationTools/NuVld/GuiDBHandler.h"
#include "ValidationTools/NuVld/GuiXmlFileHandler.h"

#ifndef _NUVLD_DB_UTILS_H_
#include "ValidationTools/NuVld/DBUtils.hh"
#endif

using std::string;
using std::ostringstream;
using std::vector;

using namespace genie;
using namespace genie::utils::str;
using namespace genie::nuvld;
using namespace genie::nuvld::constants;

ClassImp(NuVldMainFrame)

//______________________________________________________________________________
NuVldMainFrame::NuVldMainFrame(
        const TGWindow * p, UInt_t w, UInt_t h, const NuVldConfig & my_config) :
TGMainFrame(p, w, h, kFitHeight | kFitWidth)
{
  fMyConfig = new NuVldConfig(my_config);

  UInt_t kv = kVerticalFrame;
  UInt_t kh = kHorizontalFrame;
  UInt_t kf = kFitHeight | kFitWidth;

  fMain = new TGMainFrame(p,w,h,kf);
  fMain->Connect("CloseWindow()",
                         "genie::nuvld::NuVldMainFrame", this, "CloseWindow()");
  this->Init();
  this->InitializeHandlers();    // initialize GUI event handlers

  //-- define TGLayoutHints

  this->DefineLayoutHints();

  //-- add menu bar

  fMenuDock = new TGDockableFrame(fMain);
  fMenuDock->SetWindowName("NuVld Menu Bar");
  fMenuDock->EnableUndock(kTRUE);
  fMenuDock->EnableHide(kTRUE);

  fMenu = this->BuildMenuBar();

  fMenuDock->AddFrame(fMenu, fMenuBarLt);
  fMain->AddFrame(fMenuDock, new TGLayoutHints(kLHintsExpandX, 0, 0, 1, 0));

  //-- instantiate main frames (below menu & above the status bar)

  fMainFrame       = new TGCompositeFrame(fMain,       1, 1, kv|kf);
  fMainTopFrame    = new TGCompositeFrame(fMainFrame,  3, 3, kh|kf);
  fMainMiddleFrame = new TGCompositeFrame(fMainFrame,  3, 3, kh|kf);
  fMainBottomFrame = new TGCompositeFrame(fMainFrame,  3, 3, kh|kf);

  //-- TOP FRAME: add image buttons frame

  fImgBtnGrpFrm = this->BuildUpperButtonFrame();

  fMainTopFrame -> AddFrame( fImgBtnGrpFrm  );

  //-- MIDDLE FRAME:
  //         left  : neutrino / electron scattering data SQL tabs
  //         right : plotter / data viewer / session log tabs

  fMainLeftFrame   = new TGCompositeFrame(fMainMiddleFrame, 3, 3, kv);
  fMainRightFrame  = new TGCompositeFrame(fMainMiddleFrame, 3, 3, kv);

  fDBSelDock = new TGDockableFrame(fMainLeftFrame);
  fDBSelDock->SetWindowName("NuVld DB Selections");
  fDBSelDock->EnableUndock(kTRUE);
  fDBSelDock->EnableHide(kTRUE);

  fTabSql  = this->BuildSqlTab();
  fTabData = this->BuildDataTab();

  fDBSelDock        -> AddFrame ( fTabSql, fSqlTabLt   );
  fMainRightFrame   -> AddFrame ( fTabData, fDataTabLt );
  fMainLeftFrame    -> AddFrame ( fDBSelDock           );
  fMainMiddleFrame  -> AddFrame ( fMainLeftFrame,    fMLeftFrameLt  );
  fMainMiddleFrame  -> AddFrame ( fMainRightFrame,   fMRightFrameLt );

  //-- BOTTOM FRAME:
  //         left  : checkboxes
  //         right : common buttons / progress bar / selection stacking

  fBottomFrmDock = new TGDockableFrame(fMainBottomFrame);
  fBottomFrmDock->SetWindowName("NuVld");
  fBottomFrmDock->EnableUndock(kTRUE);
  fBottomFrmDock->EnableHide(kTRUE);

  fMainDBottomFrame = new TGCompositeFrame(fBottomFrmDock,    3, 3, kh|kf);
  fBottomLeftFrame  = new TGCompositeFrame(fMainDBottomFrame, 3, 3, kv);
  fBottomRightFrame = new TGCompositeFrame(fMainDBottomFrame, 3, 3, kv);

  this->AddCommonCheckButtons();                    // build second button-bar
  fStackHFrm = this->BuildSelectionStackFrame();    // add selection-stacking GUI
  fProgressBarHFrm = this->BuildLowerButtonFrame(); // add progress bar

  fBottomRightFrame ->AddFrame ( fStackHFrm,       fSelStackLt    );
  fBottomRightFrame ->AddFrame ( fProgressBarHFrm, fProgressBarLt );

  //fMainBottomFrame->AddFrame(fBottomLeftFrame,  new TGLayoutHints(kLHintsExpandX, 0, 0, 1, 0));
  fMainDBottomFrame->AddFrame(fBottomLeftFrame);
  fMainDBottomFrame->AddFrame(fBottomRightFrame, new TGLayoutHints(kLHintsExpandX, 0, 0, 1, 0));

  fBottomFrmDock->AddFrame(fMainDBottomFrame, new TGLayoutHints(kLHintsExpandX, 0, 0, 1, 0));
  fMainBottomFrame->AddFrame(fBottomFrmDock);

  //-- add top/bottom main frames to main frame, and main frame to main

  fMainFrame -> AddFrame ( fMainTopFrame    );
  fMainFrame -> AddFrame ( fMainMiddleFrame );
  fMainFrame -> AddFrame ( fMainBottomFrame );
  fMain      -> AddFrame ( fMainFrame       , new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 1, 0));

  //-- add Status Bar

  fStatusBar = this->BuildStatusBar();
  fMain->AddFrame(fStatusBar, fStatusBarLt);

  //-- initialize

  this->ResetSqlSelections();    // set all selection widgets to default values
  this->ResetFitterTab();        // set all selections to default values

  this->InitializeSyslog();      // initialize "system-log" singleton
  this->InitializeBrowser();     // initialize "data-browser" singleton
  this->ConfigHandlers();

  this->HandleDisabledNeuGEN();  // deactivate NeuGEN GUI is interface is dummy

  fMain->SetWindowName("GENIE/NuValidator");
  fMain->MapSubwindows();
  fMain->Resize( fMain->GetDefaultSize() );
  fMain->MapWindow();

  fW0 = fMain->GetWidth();  // original window width
  fH0 = fMain->GetHeight(); // original window height

  this->InitNotify();
}
//______________________________________________________________________________
NuVldMainFrame::~NuVldMainFrame()
{
  fMain->Cleanup();
  delete fMain;
  delete fMyConfig;
}
//______________________________________________________________________________

//______________________________________________________________________________
//----------------- Methods for GUI window construction ------------------------
//______________________________________________________________________________

//______________________________________________________________________________
void NuVldMainFrame::DefineLayoutHints(void)
{
  ELayoutHints kht  = kLHintsTop;
  ELayoutHints khl  = kLHintsLeft;
  ELayoutHints khr  = kLHintsRight;
  ELayoutHints khb  = kLHintsBottom;
  ELayoutHints khcx = kLHintsCenterX;
  ELayoutHints khcy = kLHintsCenterY;
  ELayoutHints khex = kLHintsExpandX;
  ELayoutHints khey = kLHintsExpandY;

  ULong_t hintMenuBarLayout       = kht  | khl  | khex;
  ULong_t hintMenuBarItemLayout   = kht  | khl;
  ULong_t hintMenuBarHelpLayout   = kht  | khr;
  ULong_t hintTabPlotterLayout    = kht  | khl  | khex | khey;
  ULong_t hintTabFitterLayout     = kht  | khl  | khex | khey;
  ULong_t hintTabLogLayout        = kht  | khl  | khex | khey;
  ULong_t hintTabDataViewerLayout = kht  | khl  | khex | khey;
  ULong_t hintTabNuSqlLayout      = kht  | khl  | khex | khey;
  ULong_t hintTabElSqlLayout      = kht  | khl  | khex | khey;
  ULong_t hintTabSFSqlLayout      = kht  | khl  | khex | khey;
  ULong_t hintTabDataLayout       = kht  | khex | khey;
  ULong_t hintTabSqlLayout        = kht  | khex | khey;
  ULong_t hintExitButtonLayout    = khr;
  ULong_t hintLeftButtonsLayout   = khl;
  ULong_t hintStatusBarLayout     = khb  | khl  | khex;
  ULong_t hintMLeftFrameLayout    = khcy;
  ULong_t hintMRightFrameLayout   = kht  | khex | khey;
  ULong_t hintProgressBarLayout   = khcx;
  ULong_t hintSelectStackLayout   = khcx;
  ULong_t hintFitLeftFrameLayout  = kht  | khl;
  ULong_t hintFitRightFrameLayout = kht  | khr;

  fMenuBarLt       = new TGLayoutHints(hintMenuBarLayout,       0, 0,  1, 1);
  fMenuBarItemLt   = new TGLayoutHints(hintMenuBarItemLayout,   0, 4,  0, 0);
  fMenuBarHelpLt   = new TGLayoutHints(hintMenuBarHelpLayout);
  fPlotterTabLt    = new TGLayoutHints(hintTabPlotterLayout,    5, 5, 10, 1);
  fLogTabLt        = new TGLayoutHints(hintTabLogLayout,        5, 5, 10, 1);
  fFitterTabLt     = new TGLayoutHints(hintTabFitterLayout,     5, 5, 10, 1);
  fDataViewTabLt   = new TGLayoutHints(hintTabDataViewerLayout, 5, 5, 10, 1);
  fNuSqlTabLt      = new TGLayoutHints(hintTabNuSqlLayout,      5, 5, 10, 1);
  fElSqlTabLt      = new TGLayoutHints(hintTabElSqlLayout,      5, 5, 10, 1);
  fSFSqlTabLt      = new TGLayoutHints(hintTabSFSqlLayout,      5, 5, 10, 1);
  fDataTabLt       = new TGLayoutHints(hintTabDataLayout,       5, 5, 10, 1);
  fSqlTabLt        = new TGLayoutHints(hintTabSqlLayout,        5, 5, 10, 1);
  fProgressBarLt   = new TGLayoutHints(hintProgressBarLayout,   5, 5,  3, 4);
  fExitBtnLt       = new TGLayoutHints(hintExitButtonLayout,    5, 5,  3, 4);
  fLeftBtnLt       = new TGLayoutHints(hintLeftButtonsLayout,   5, 5,  3, 4);
  fStatusBarLt     = new TGLayoutHints(hintStatusBarLayout,     0, 0,  2, 0);
  fMLeftFrameLt    = new TGLayoutHints(hintMLeftFrameLayout,    1, 1,  1, 1);
  fMRightFrameLt   = new TGLayoutHints(hintMRightFrameLayout,   1, 1,  1, 1);
  fProgressBarLt   = new TGLayoutHints(hintProgressBarLayout,   2, 2,  2, 2);
  fSelStackLt      = new TGLayoutHints(hintSelectStackLayout,   2, 2,  2, 2);
  fFitLeftFrameLt  = new TGLayoutHints(hintFitLeftFrameLayout,  1, 1,  1, 1);
  fFitRightFrameLt = new TGLayoutHints(hintFitRightFrameLayout, 1, 1,  1, 1);
}
//______________________________________________________________________________
TGMenuBar * NuVldMainFrame::BuildMenuBar(void)
{
  // Menu: File

  fMenuFile = new TGPopupMenu(gClient->GetRoot());

  fMenuFile->AddEntry ("&Open XML",               M_FILE_OPEN  );
  fMenuFile->AddEntry ("&Close XML",              M_FILE_CLOSE );
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry ("&Test XML File Parsing",  M_FILE_PARSE );
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry ("&Load Selections",        M_LOAD_STACKED_DATA );
  fMenuFile->AddEntry ("&Save Selections",        M_SAVE_STACKED_DATA );
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry ("E&xit",                   M_FILE_EXIT  );

  fMenuFile->DisableEntry(M_FILE_CLOSE); // disable "close" until next "open"

  fMenuFile->Connect("Activated(Int_t)",
                                  "genie::nuvld::NuVldMainFrame",
                                                     this,"HandleMenu(Int_t)");

  fXmlFileHandler->AttachFileMenu(fMenuFile); // to (de/)activate open/close

    // Menu: View

  fMenuView = new TGPopupMenu(gClient->GetRoot());

  fMenuView->AddEntry ("Enable (Un)Dock Menu-Bar",       M_VIEW_ENABLE_UNDOCK_MENU           );
  fMenuView->AddEntry ("Enable (Un)Dock Selection-Tabs", M_VIEW_ENABLE_UNDOCK_SELECTION_TABS );
  fMenuView->AddEntry ("Enable (Un)Dock Bottom-Frame",   M_VIEW_ENABLE_UNDOCK_BOTTOM_FRAME   );
  fMenuView->AddEntry ("Enable Hide Menu-Bar",           M_VIEW_ENABLE_HIDE_MENU             );
  fMenuView->AddEntry ("Enable Hide Selection-Tabs",     M_VIEW_ENABLE_HIDE_SELECTION_TABS   );
  fMenuView->AddEntry ("Enable Hide Bottom-Frame",       M_VIEW_ENABLE_HIDE_BOTTOM_FRAME     );
  fMenuView->AddSeparator();
  fMenuView->AddEntry ("Dock Menu-Bar",                  M_VIEW_DOCK_MENU                    );
  fMenuView->AddEntry ("UnDock Menu-Bar",                M_VIEW_UNDOCK_MENU                  );
  fMenuView->AddEntry ("Hide Menu-Bar",                  M_VIEW_HIDE_MENU                    );
  fMenuView->AddSeparator();
  fMenuView->AddEntry ("Dock Selection-Tabs",            M_VIEW_DOCK_SELECTION_TABS          );
  fMenuView->AddEntry ("UnDock Selection-Tabs",          M_VIEW_UNDOCK_SELECTION_TABS        );
  fMenuView->AddEntry ("Hide Selection-Tabs",            M_VIEW_HIDE_SELECTION_TABS          );
  fMenuView->AddEntry ("Show Selection-Tabs",            M_VIEW_SHOW_SELECTION_TABS          );
  fMenuView->AddSeparator();
  fMenuView->AddEntry ("Dock Bottom-Frame",              M_VIEW_DOCK_BOTTOM_FRAME            );
  fMenuView->AddEntry ("UnDock Bottom-Frame",            M_VIEW_UNDOCK_BOTTOM_FRAME          );
  fMenuView->AddEntry ("Hide Bottom-Frame",              M_VIEW_HIDE_BOTTOM_FRAME            );
  fMenuView->AddEntry ("Show Bottom-Frame",              M_VIEW_SHOW_BOTTOM_FRAME            );

  fMenuView->CheckEntry ( M_VIEW_ENABLE_UNDOCK_MENU           );
  fMenuView->CheckEntry ( M_VIEW_ENABLE_HIDE_MENU             );
  fMenuView->CheckEntry ( M_VIEW_ENABLE_UNDOCK_SELECTION_TABS );
  fMenuView->CheckEntry ( M_VIEW_ENABLE_HIDE_SELECTION_TABS   );
  fMenuView->CheckEntry ( M_VIEW_ENABLE_UNDOCK_BOTTOM_FRAME   );
  fMenuView->CheckEntry ( M_VIEW_ENABLE_HIDE_BOTTOM_FRAME     );

  fMenuView->Connect("Activated(Int_t)",
                                  "genie::nuvld::NuVldMainFrame",
                                                     this,"HandleMenu(Int_t)");
  // Menu: DataBase

  fMenuDBase = new TGPopupMenu(gClient->GetRoot());

  fMenuDBase->AddEntry ("&Connect",             M_DATA_CONNECT    );
  fMenuDBase->AddEntry ("Disconnect",           M_DATA_CLOSE      );
  fMenuDBase->AddEntry ("Check connection",     M_DATA_VERIFY     );
  fMenuDBase->AddEntry ("DBase info",           M_DATA_DBINFO     );
  fMenuDBase->AddSeparator();
  fMenuDBase->AddEntry ("&Bootstrap",           M_DATA_BOOTSTRAP  );
  fMenuDBase->AddEntry ("&UpLoad XML",          M_DATA_UPLOAD     );
  fMenuDBase->AddSeparator();
  fMenuDBase->AddEntry ("Type SQL query",       M_DATA_QUERY      );
  fMenuDBase->AddEntry ("Load SQL query",       M_DATA_QUERY_FILE );
  fMenuDBase->AddSeparator();
  fMenuDBase->AddEntry ("Plot data",       M_DATA_QUERY_DRAW_GUI  );
  fMenuDBase->AddEntry ("Print data",     M_DATA_QUERY_PRINT_GUI  );

  fMenuDBase->Connect("Activated(Int_t)",
                       "genie::nuvld::NuVldMainFrame",this,"HandleMenu(Int_t)");

  // Menu: Export/Import

  fMenuExpImp = new TGPopupMenu(gClient->GetRoot());

  fMenuExpImp->AddEntry("Plot  -> eps/gif plot or ROOT macro", M_EXPORT_PLOT);
  fMenuExpImp->AddEntry("Table -> text file",                  M_EXPORT_TABLE);
  fMenuExpImp->AddEntry("Model -> XML/ROOT/text file",         M_EXPORT_MODEL);
  fMenuExpImp->AddSeparator();
  fMenuExpImp->AddEntry("XML/text file -> Model",              M_IMPORT_MODEL);

  fMenuExpImp->Connect("Activated(Int_t)",
                                   "genie::nuvld::NuVldMainFrame",
                                                    this,"HandleMenu(Int_t)");
/*
  // Menu: NeuGEN

  fMenuNeuGen = new TGPopupMenu(gClient->GetRoot());

  fMenuNeuGen->AddEntry("Physics",       M_NEUGEN_CONFIG_PHYSICS);
  fMenuNeuGen->AddEntry("Process",       M_NEUGEN_CONFIG_PROCESS);
  fMenuNeuGen->AddEntry("Run",           M_NEUGEN_RUN);

  fMenuNeuGen->Connect("Activated(Int_t)",
                                   "genie::nuvld::NuVldMainFrame",
                                                 this,"HandleMenu(Int_t)");
*/

  // Menu: GENIE

  fMenuGENIE = new TGPopupMenu(gClient->GetRoot());

  // Menu: Fit

  fMenuFit = new TGPopupMenu(gClient->GetRoot());

  fMenuFit->AddEntry("Open",  M_FIT_OPEN );
  fMenuFit->AddEntry("Run",   M_FIT_RUN  );
  fMenuFit->AddEntry("Reset", M_FIT_RESET);
  fMenuFit->Connect("Activated(Int_t)",
                                "genie::nuvld::NuVldMainFrame",
                                                    this,"HandleMenu(Int_t)");
  // Menu: Help

  fMenuHelp = new TGPopupMenu(gClient->GetRoot());

  fMenuHelp->AddEntry("&About",           M_HELP_ABOUT );
  fMenuHelp->AddEntry("www:NuValidator",  M_HELP_WWW   );
  fMenuHelp->AddEntry("www::DURHAM",      M_HELP_DURHAM);
  fMenuHelp->AddSeparator();
  fMenuHelp->AddEntry("Connecting to the dbase", M_HELP_HOWTO_CONN_DBASE);
  fMenuHelp->AddEntry("How to fill the dbase",   M_HELP_HOWTO_FILL_DBASE);

  fMenuHelp->Connect("Activated(Int_t)",
                                "genie::nuvld::NuVldMainFrame",
                                                    this,"HandleMenu(Int_t)");

  //menu_bar = new TGMenuBar(fMain, 1, 1, kHorizontalFrame);
  TGMenuBar * menu_bar = new TGMenuBar(fMenuDock, 1, 1, kHorizontalFrame);

  menu_bar->AddPopup("&File",          fMenuFile,   fMenuBarItemLt);
  menu_bar->AddPopup("&View",          fMenuView,   fMenuBarItemLt);
  menu_bar->AddPopup("&Database",      fMenuDBase,  fMenuBarItemLt);
  menu_bar->AddPopup("&Export/Import", fMenuExpImp, fMenuBarItemLt);
//menu_bar->AddPopup("&NeuGEN",        fMenuNeuGen, fMenuBarItemLt);
  menu_bar->AddPopup("&GENIE",         fMenuGENIE,  fMenuBarItemLt);
  menu_bar->AddPopup("&Fit",           fMenuFit,    fMenuBarItemLt);
  menu_bar->AddPopup("&Help",          fMenuHelp,   fMenuBarHelpLt);

  return menu_bar;
}
//______________________________________________________________________________
TGGroupFrame * NuVldMainFrame::BuildUpperButtonFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(fMainTopFrame, "", kVerticalFrame);

  fBtnMatrixLt = new TGMatrixLayout(grpf, 0, 15, 1);

  grpf->SetLayoutManager( fBtnMatrixLt );

  this->CreateUpperFrameButtons(grpf);   // create buttons
  this->SetUpperFrameButtonText();       // set tool-tip text
  this->ConnectUpperFrameButtons();      // connect Click() with a method

  grpf -> AddFrame( fOpenXmlBtn      );
  grpf -> AddFrame( fParseXmlBtn     );
  grpf -> AddFrame( fDBConnectBtn    );
  grpf -> AddFrame( fDBCloseBtn      );
  grpf -> AddFrame( fDBInfoBtn       );
  grpf -> AddFrame( fDBCheckBtn      );
  grpf -> AddFrame( fDBBootstrapBtn  );
  grpf -> AddFrame( fSqlQInpBtn      );
  grpf -> AddFrame( fSqlQFileBtn     );
  grpf -> AddFrame( fDBUploadBtn     );
  grpf -> AddFrame( fNeugenConfigBtn );
  grpf -> AddFrame( fNeugenProcBtn   );
  grpf -> AddFrame( fHelpBtn         );
  grpf -> AddFrame( fDurhamBtn       );
  grpf -> AddFrame( fAboutBtn        );

  return grpf;
}
//______________________________________________________________________________
void NuVldMainFrame::CreateUpperFrameButtons(TGGroupFrame * gf)
{
  fOpenXmlBtn      = new TGPictureButton(gf, this->Pic("open",         32,32));
  fParseXmlBtn     = new TGPictureButton(gf, this->Pic("parse",        32,32));
  fDBConnectBtn    = new TGPictureButton(gf, this->Pic("connect",      32,32));
  fDBCloseBtn      = new TGPictureButton(gf, this->Pic("dbclose",      32,32));
  fDBInfoBtn       = new TGPictureButton(gf, this->Pic("print_dbinfo", 32,32));
  fDBCheckBtn      = new TGPictureButton(gf, this->Pic("dbcheck",      32,32));
  fDBBootstrapBtn  = new TGPictureButton(gf, this->Pic("bootstrap",    32,32));
  fSqlQInpBtn      = new TGPictureButton(gf, this->Pic("sql",          32,32));
  fSqlQFileBtn     = new TGPictureButton(gf, this->Pic("sql_file",     32,32));
  fDBUploadBtn     = new TGPictureButton(gf, this->Pic("upload",       32,32));
  fNeugenConfigBtn = new TGPictureButton(gf, this->Pic("neugen_config",32,32));
  fNeugenProcBtn   = new TGPictureButton(gf, this->Pic("neugen_inputs",32,32));
  fHelpBtn         = new TGPictureButton(gf, this->Pic("nuvld_help",   32,32));
  fDurhamBtn       = new TGPictureButton(gf, this->Pic("durham",       32,32));
  fAboutBtn        = new TGPictureButton(gf, this->Pic("about",        32,32));
}
//______________________________________________________________________________
void NuVldMainFrame::SetUpperFrameButtonText(void)
{
  fOpenXmlBtn      -> SetToolTipText( "Open XML file" ,                1);
  fParseXmlBtn     -> SetToolTipText( "Parse XML file",                1);
  fDBConnectBtn    -> SetToolTipText( "Connect to dbase",              1);
  fDBCloseBtn      -> SetToolTipText( "Disconnect from dbase",         1);
  fDBInfoBtn       -> SetToolTipText( "Print dbase info",              1);
  fDBCheckBtn      -> SetToolTipText( "Check dbase connection",        1);
  fDBBootstrapBtn  -> SetToolTipText( "Bootsrap dbase",                1);
  fSqlQInpBtn      -> SetToolTipText( "Send SQL Query",                1);
  fSqlQFileBtn     -> SetToolTipText( "Send custom SQL Query",         1);
  fDBUploadBtn     -> SetToolTipText( "Upload XML to dbase",           1);
  fNeugenConfigBtn -> SetToolTipText( "Config. NeuGEN physics input",  1);
  fNeugenProcBtn   -> SetToolTipText( "Config. NeuGEN input card",     1);
  fHelpBtn         -> SetToolTipText( "NuValidator Online Help",       1);
  fDurhamBtn       -> SetToolTipText( "DURHAM Neutrino Web Page",      1);
  fAboutBtn        -> SetToolTipText( "About NuValidator",             1);
}
//______________________________________________________________________________
void NuVldMainFrame::ConnectUpperFrameButtons(void)
{
  fOpenXmlBtn      -> Connect("Clicked()","genie::nuvld::GuiXmlFileHandler",
                                                fXmlFileHandler,"OpenFile()");
  fParseXmlBtn     -> Connect("Clicked()","genie::nuvld::GuiXmlFileHandler",
                                               fXmlFileHandler,"ParseFile()");
  fDBConnectBtn    -> Connect("Clicked()","genie::nuvld::GuiDBHandler",
                                           fDBaseHandler, "MakeConnection()");
  fDBCloseBtn      -> Connect("Clicked()","genie::nuvld::GuiDBHandler",
                                           fDBaseHandler,"CloseConnection()");
  fDBInfoBtn       -> Connect("Clicked()","genie::nuvld::GuiDBHandler",
                                                 fDBaseHandler,"PrintInfo()");
  fDBCheckBtn      -> Connect("Clicked()","genie::nuvld::GuiDBHandler",
                                           fDBaseHandler,"CheckConnection()");
  fDBBootstrapBtn  -> Connect("Clicked()","genie::nuvld::GuiDBHandler",
                                                 fDBaseHandler,"Bootstrap()");
  fSqlQInpBtn      -> Connect("Clicked()","genie::nuvld::GuiDBHandler",
                                    fDBaseHandler,"QueryWithSqlFromDialog()");
  fSqlQFileBtn     -> Connect("Clicked()","genie::nuvld::GuiDBHandler",
                                      fDBaseHandler,"QueryWithSqlFromFile()");
  fDBUploadBtn     -> Connect("Clicked()","genie::nuvld::GuiXmlFileHandler",
                                             fXmlFileHandler,"SendToDBase()");
  fNeugenProcBtn   -> Connect("Clicked()","genie::nuvld::NuVldMainFrame",
                                                this,"ConfigNeugenProcess()");
  fNeugenConfigBtn -> Connect("Clicked()","genie::nuvld::NuVldMainFrame",
                                                this,"ConfigNeugenPhysics()");
  fHelpBtn         -> Connect("Clicked()","genie::nuvld::GuiHelpHandler",
                                                fHelpHandler,"NuVldOnline()");
  fDurhamBtn       -> Connect("Clicked()","genie::nuvld::GuiHelpHandler",
                                               fHelpHandler,"DurhamOnline()");
  fAboutBtn        -> Connect("Clicked()","genie::nuvld::GuiHelpHandler",
                                                 fHelpHandler,"NuVldAbout()");
}
//______________________________________________________________________________
TGTab * NuVldMainFrame::BuildSqlTab(void)
{
  TGCompositeFrame * tf = 0;

  unsigned int width  = 250;
  unsigned int height = 550;

  //TGTab * tab = new TGTab(fMainLeftFrame, 1, 1);
  TGTab * tab = new TGTab(fDBSelDock, 1, 1);

  //-- tab: SQL GUI widgets for v scattering data

  tf = tab->AddTab( "vN" );

  fTabNuSql = fNuXSecTab->Create(tf, width, height);
  tf -> AddFrame( fTabNuSql, fNuSqlTabLt );

  //-- tab: SQL GUI widgets for e scattering data

  tf = tab->AddTab( "eN" );

  fTabElSql = fElXSecTab->Create(tf, width, height);
  tf   -> AddFrame( fTabElSql, fElSqlTabLt );

  //-- tab: SQL GUI widgets for e scattering data

  tf = tab->AddTab( "S/F" );

  fTabSFSql = fSFTab->Create(tf, width, height);
  tf -> AddFrame( fTabSFSql, fSFSqlTabLt );

  return tab;
}
//______________________________________________________________________________
void NuVldMainFrame::AddCommonCheckButtons(void)
{
  TGCompositeFrame * frm = fBottomLeftFrame;

  fShowColorCodeChkB = new TGCheckButton(frm, "Color-Code plot",  72);
  fShowExtLegendChkB = new TGCheckButton(frm, "External legend",  73);
  fUseStackedChkB    = new TGCheckButton(frm, "Use stacked",      74);

  frm -> AddFrame( fShowColorCodeChkB );
  frm -> AddFrame( fShowExtLegendChkB );
  frm -> AddFrame( fUseStackedChkB    );
}
//______________________________________________________________________________
TGTab * NuVldMainFrame::BuildDataTab(void)
{
// Build tabs

  TGCompositeFrame * tf = 0;

  unsigned int width  = 650;
  unsigned int height = 550;

  TGTab * tab = new TGTab(fMainRightFrame, 1, 1);
  //TGTab * tab = new TGTab(fMainTabDock, 1, 1);

  //--- tab: "Plotter"

  tf = tab->AddTab( "Plotter" );

  fTabPlotter    = new TGCompositeFrame    (tf, width, height, kVerticalFrame);
  fPlotTabEmbCnv = new TRootEmbeddedCanvas ("evanvas", fTabPlotter, width, height);

  fPlotTabEmbCnv -> GetCanvas() -> SetBorderMode (0);
  fPlotTabEmbCnv -> GetCanvas() -> SetFillColor  (0);

  fTabPlotter -> AddFrame( fPlotTabEmbCnv,    fPlotterTabLt );
  tf          -> AddFrame( fTabPlotter, fPlotterTabLt );

  //--- tab "Data Viewer"

  tf = tab->AddTab("Data Viewer");

  fTabDataViewer = new TGCompositeFrame(tf,  width, height, kVerticalFrame);

  fDataViewer = new TGTextEdit(fTabDataViewer,  width, height, kSunkenFrame | kDoubleBorder);
  fDataViewer->AddLine( "Data-Base Viewer:" );

  fTabDataViewer -> AddFrame( fDataViewer,    fDataViewTabLt);
  tf             -> AddFrame( fTabDataViewer, fDataViewTabLt);

  //--- tab: "Fitter"

  tf = tab->AddTab( "Fitter" );

  fTabFitter  = new TGCompositeFrame(tf, width, height, kHorizontalFrame);

  this->FillFitterFrame();

  tf->AddFrame( fTabFitter, fFitterTabLt );

  //--- tab: "Session Log"

  tf = tab->AddTab("Session Log");

  fTabLog = new TGCompositeFrame(tf, width, height, kVerticalFrame);

  fLog = new TGTextEdit(fTabLog,  width, height, kSunkenFrame | kDoubleBorder);
  fLog->AddLine( "NuValidator Log:" );

  fTabLog -> AddFrame( fLog,    fLogTabLt);
  tf      -> AddFrame( fTabLog, fLogTabLt);

  return tab;
}
//______________________________________________________________________________
void NuVldMainFrame::FillFitterFrame(void)
{
  // left + right fitter frames

  fFitterLeftFrame   = new TGCompositeFrame(fTabFitter, 3, 3, kVerticalFrame);
  fFitterRightFrame  = new TGCompositeFrame(fTabFitter, 3, 3, kVerticalFrame);

  fLFitSpacerLb  = new TGLabel( fFitterLeftFrame,  new TGString(" "));
  fRFitSpacerLb  = new TGLabel( fFitterRightFrame, new TGString(" "));

  // fitter frame: fitter combo

  fFitterGrpFrm = new TGGroupFrame(fFitterLeftFrame,
                                    "Select Fitter/Params", kHorizontalFrame);

  fFitterCBx = new TGComboBox(fFitterGrpFrm, 112);

  utils::gui::FillComboBox( fFitterCBx, kFitters );

  fFitterCBx->Resize (160,20);

  fSelectNeuGenFitParams = new TGTextButton (
                                  fFitterGrpFrm, "Select fit parameters", 177);

  fSelectNeuGenFitParams->Connect("Clicked()",
              "genie::nuvld::NuVldMainFrame", this, "SelectNeuGenFitParams()");

  // add option to set a range for the free parameters (E, W, Q^2...)

  fFitFreeParamGrpFrm = new TGGroupFrame(
                     fFitterLeftFrame,"Free parameter: range", kVerticalFrame);

  fFitFreeParamGrpFrm->SetTitlePos(TGGroupFrame::kLeft);
  fFitFreeParamGrpFrm->SetLayoutManager(
                             new TGMatrixLayout(fFitFreeParamGrpFrm, 0, 4, 1));

  fXMinNmE = new TGNumberEntry(
                       fFitFreeParamGrpFrm, 0, 6, 1, TGNumberFormat::kNESReal);
  fXMaxNmE = new TGNumberEntry(
                     fFitFreeParamGrpFrm, 100, 6, 1, TGNumberFormat::kNESReal);

  fXMinLb = new TGLabel(fFitFreeParamGrpFrm, new TGString( " E-min: "));
  fXMaxLb = new TGLabel(fFitFreeParamGrpFrm, new TGString( " E-max: "));

  // fit: picture buttons

  fFitBtnGrpFrm = new TGGroupFrame(
             fFitterRightFrame,"Fitter: shortcut buttons", kHorizontalFrame);

  fDoFitBtn     = new TGPictureButton(fFitBtnGrpFrm, this->Pic("fit",    32,32));
  fPrmScanBtn   = new TGPictureButton(fFitBtnGrpFrm, this->Pic("scan",   32,32));
  fPrmScan1dBtn = new TGPictureButton(fFitBtnGrpFrm, this->Pic("scan1d", 32,32));
  fPrmScan2dBtn = new TGPictureButton(fFitBtnGrpFrm, this->Pic("scan2d", 32,32));
  fResetFitBtn  = new TGPictureButton(fFitBtnGrpFrm, this->Pic("reset",  32,32));

  fDoFitBtn     -> SetToolTipText( "Do the fit",                              1);
  fPrmScanBtn   -> SetToolTipText( "multi-D parameter space MC scanning",     2);
  fPrmScan1dBtn -> SetToolTipText( "1-D parameter space scanning",            3);
  fPrmScan2dBtn -> SetToolTipText( "2-D parameter space scanning",            4);
  fResetFitBtn  -> SetToolTipText( "Reset fitter selections" ,                5);

  fDoFitBtn    -> Connect("Clicked()","genie::nuvld::NuVldMainFrame",
                                                          this,"RunFitter()");
  fPrmScanBtn  -> Connect("Clicked()","genie::nuvld::NuVldMainFrame",
                                                       this,"RunMcScanner()");
  fPrmScan1dBtn-> Connect("Clicked()","genie::nuvld::NuVldMainFrame",
                                                       this,"Run1dScanner()");
  fPrmScan2dBtn-> Connect("Clicked()","genie::nuvld::NuVldMainFrame",
                                                       this,"Run2dScanner()");
  fResetFitBtn -> Connect("Clicked()", "genie::nuvld::NuVldMainFrame",
                                                     this,"ResetFitterTab()");

  // text results TGTextEdit

  unsigned int width  = 330;
  unsigned int height = 230;

  fFitTxtResults = new TGTextEdit(fFitterLeftFrame, width, 2*height,
                                                kSunkenFrame | kDoubleBorder);

  fFitTxtResults->AddLine( "Fit Results:" );

  // fitter pad embedded canvasses

  //- data / fitted func

  fFitTabFuncEmbCnv = new TRootEmbeddedCanvas ("evanvas_small",
                                          fFitterRightFrame, width, height);

  fFitTabFuncEmbCnv -> GetCanvas() -> SetBorderMode (0);
  fFitTabFuncEmbCnv -> GetCanvas() -> SetFillColor  (0);

  //- chisq

  fFitTabChisqEmbCnv = new TRootEmbeddedCanvas ("evanvas_chisq",
                                          fFitterRightFrame, width, height);

  fFitTabChisqEmbCnv -> GetCanvas() -> SetBorderMode (0);
  fFitTabChisqEmbCnv -> GetCanvas() -> SetFillColor  (0);

  // arrange frames

  fFitterGrpFrm -> AddFrame( fFitterCBx );
  fFitterGrpFrm -> AddFrame( fSelectNeuGenFitParams );

  fFitFreeParamGrpFrm -> AddFrame( fXMinLb  );
  fFitFreeParamGrpFrm -> AddFrame( fXMinNmE );
  fFitFreeParamGrpFrm -> AddFrame( fXMaxLb  );
  fFitFreeParamGrpFrm -> AddFrame( fXMaxNmE );

  fFitBtnGrpFrm -> AddFrame( fDoFitBtn     );
  fFitBtnGrpFrm -> AddFrame( fPrmScanBtn   );
  fFitBtnGrpFrm -> AddFrame( fPrmScan1dBtn );
  fFitBtnGrpFrm -> AddFrame( fPrmScan2dBtn );
  fFitBtnGrpFrm -> AddFrame( fResetFitBtn  );

  fFitterLeftFrame -> AddFrame( fFitterGrpFrm       );
  fFitterLeftFrame -> AddFrame( fFitFreeParamGrpFrm );
  fFitterLeftFrame -> AddFrame( fLFitSpacerLb       );
  fFitterLeftFrame -> AddFrame( fFitTxtResults      );

  fFitterRightFrame -> AddFrame( fFitBtnGrpFrm      );
  fFitterRightFrame -> AddFrame( fRFitSpacerLb      );
  fFitterRightFrame -> AddFrame( fFitTabFuncEmbCnv  );
  fFitterRightFrame -> AddFrame( fRFitSpacerLb      );
  fFitterRightFrame -> AddFrame( fFitTabChisqEmbCnv );

  fTabFitter -> AddFrame ( fFitterLeftFrame,  fFitLeftFrameLt  );
  fTabFitter -> AddFrame ( fFitterRightFrame, fFitRightFrameLt );
}
//______________________________________________________________________________
void NuVldMainFrame::SelectNeuGenFitParams(void)
{
/*
  new NeuGenFitParamsDialog(
                    gClient->GetRoot(), fMain, 380, 250, kVerticalFrame, fNGFP);
*/
}
//______________________________________________________________________________
TGHorizontalFrame * NuVldMainFrame::BuildSelectionStackFrame(void)
{
  TGHorizontalFrame * hf = new TGHorizontalFrame(fBottomRightFrame, 10, 10);

  // selection-stack-frame: labels

  fStackDBTableLb = new TGLabel(hf, new TGString("STACK: Data> "));
  fStackConfigLb  = new TGLabel(hf, new TGString(" Config>"));
  fLinkSelLb      = new TGLabel(hf, new TGString(" Link> "));

  // selection-stack-frame: text entries

  fStackTableNameTxE  = new TGTextEntry(hf, new TGTextBuffer(10) );
  fStackConfigNameTxE = new TGTextEntry(hf, new TGTextBuffer(10) );

  // selection-stack-frame: picture buttons

  fStackTableBtn  = new TGPictureButton(hf, this->Pic("add_stack_item",   24,24));
  fStackConfigBtn = new TGPictureButton(hf, this->Pic("add_stack_item",   24,24));
  fLinkStackedBtn = new TGPictureButton(hf, this->Pic("link_selections",  24,24));
  fDelStackedBtn  = new TGPictureButton(hf, this->Pic("delete_stack_item",24,24));

  fStackTableBtn -> SetToolTipText( "Stack a data selection", 1);
  fDelStackedBtn -> SetToolTipText( "Erase a data selection", 2);


  fStackTableBtn  -> Connect("Clicked()",
                 "genie::nuvld::GuiStackHandler",fStackHandler,"StackDBTable()");
  fStackConfigBtn -> Connect("Clicked()",
                  "genie::nuvld::GuiStackHandler",fStackHandler,"StackConfig()");
  fDelStackedBtn  -> Connect("Clicked()",
             "genie::nuvld::GuiStackHandler",fStackHandler,"EraseStackedItem()");

  // selection-stack-frame: combo boxes

  fTableStackCBx  = new TGComboBox(hf, 111);
  fConfigStackCBx = new TGComboBox(hf, 112);

  fTableStackCBx  -> Resize (100,20);
  fConfigStackCBx -> Resize (100,20);

  // selection-stack-frame: add all widgets to TGHorizontalFrame

  hf -> AddFrame ( fStackDBTableLb     );
  hf -> AddFrame ( fStackTableNameTxE  );
  hf -> AddFrame ( fStackTableBtn      );
  hf -> AddFrame ( fStackConfigLb      );
  hf -> AddFrame ( fStackConfigNameTxE );
  hf -> AddFrame ( fStackConfigBtn     );
  hf -> AddFrame ( fLinkSelLb          );
  hf -> AddFrame ( fTableStackCBx      );
  hf -> AddFrame ( fConfigStackCBx     );
  hf -> AddFrame ( fLinkStackedBtn     );
  hf -> AddFrame ( fDelStackedBtn      );

  return hf;
}
//______________________________________________________________________________
TGHorizontalFrame * NuVldMainFrame::BuildLowerButtonFrame(void)
{
  TGHorizontalFrame * hf = new TGHorizontalFrame(fBottomRightFrame, 10, 10);

  // lower frame: left hand-side picture buttons

  fSelResetBtn  = new TGPictureButton(hf, this->Pic("reset",      32,32));
  fViewClearBtn = new TGPictureButton(hf, this->Pic("viewclear",  32,32));
  fDrawDataBtn  = new TGPictureButton(hf, this->Pic("draw_data",  32,32));
  fPrintDataBtn = new TGPictureButton(hf, this->Pic("print_data", 32,32));
  fNeugenRunBtn = new TGPictureButton(hf, this->Pic("neugen_run", 32,32));
  fSaveBtn      = new TGPictureButton(hf, this->Pic("save",       32,32));

  string sfSelResetBtn  = "Reset selection";
  string sfViewClearBtn = "Clear view";
  string sfDrawDataBtn  = "Query dbase & draw data";
  string sfPrintDataBtn = "Query dbase & print data in text format";
  string sfNeugenRunBtn = "Run NeuGEN with selected inputs & draw its prediction";
  string sfSaveBtn      = "Save graph";

  fSelResetBtn  -> SetToolTipText( sfSelResetBtn.c_str(),  1 );
  fViewClearBtn -> SetToolTipText( sfViewClearBtn.c_str(), 1 );
  fDrawDataBtn  -> SetToolTipText( sfDrawDataBtn.c_str(),  1 );
  fPrintDataBtn -> SetToolTipText( sfPrintDataBtn.c_str(), 1 );
  fNeugenRunBtn -> SetToolTipText( sfNeugenRunBtn.c_str(), 1 );
  fSaveBtn      -> SetToolTipText( sfSaveBtn.c_str(),      1 );

  string tc = "genie::nuvld::NuVldMainFrame"; // this class

  fSelResetBtn  -> Connect("Clicked()",tc.c_str(),this,"ResetSqlSelections()" );
  fViewClearBtn -> Connect("Clicked()",tc.c_str(),this,"ClearViewer()"        );
  fDrawDataBtn  -> Connect("Clicked()",tc.c_str(),this,"DrawDBTable()"        );
  fPrintDataBtn -> Connect("Clicked()",tc.c_str(),this,"PrintDBTable()"       );
  fNeugenRunBtn -> Connect("Clicked()",tc.c_str(),this,"RunNeuGen()"          );
  fSaveBtn      -> Connect("Clicked()",tc.c_str(),this,"HandleSaveCanvas()"   );

  // lower frame: progress bar

  fProgressBar = new TGHProgressBar(hf, TGProgressBar::kStandard, 300);
  fProgressBar->SetRange(0., 100.);
  fProgressBar->SetBarColor("purple");
  fProgressBar->ShowPosition(true, true, "progress");
  fProgressBar->SetFillType(TGProgressBar::kBlockFill);

  // lower frame: exit button

  fExitBtn = new TGPictureButton(hf, this->Pic("exit",32,32),
                                                 "gApplication->Terminate(0)");

  fExitBtn -> SetToolTipText( "Exit GUI", 1);

  // add widgets to TGHorizontalFrame

  hf -> AddFrame ( fSelResetBtn,   fLeftBtnLt     );
  hf -> AddFrame ( fViewClearBtn,  fLeftBtnLt     );
  hf -> AddFrame ( fDrawDataBtn,   fLeftBtnLt     );
  hf -> AddFrame ( fPrintDataBtn,  fLeftBtnLt     );
  hf -> AddFrame ( fNeugenRunBtn,  fLeftBtnLt     );
  hf -> AddFrame ( fSaveBtn,       fLeftBtnLt     );
  hf -> AddFrame ( fProgressBar,   fProgressBarLt );
  hf -> AddFrame ( fExitBtn,       fExitBtnLt     );

  return hf;
}
//______________________________________________________________________________
TGStatusBar * NuVldMainFrame::BuildStatusBar(void)
{
//! Adds a status bar at the main frame

  Int_t parts[] = { 60, 20, 20 };

  TGStatusBar * status_bar = new TGStatusBar(fMain, 50, 10, kHorizontalFrame);

  status_bar->SetParts(parts, 3);

  return status_bar;
}
//______________________________________________________________________________
const char * NuVldMainFrame::Icon(const char * name)
{
//! Returns the full path-name for the requested picture

  ostringstream pic;
  pic  << gSystem->Getenv("GENIE") << "/data/icons/" << name << ".xpm";

  return pic.str().c_str();
}
//______________________________________________________________________________
const TGPicture * NuVldMainFrame::Pic(const char * name, int x, int y)
{
//! Returns the requested picture at the requested size

  return gClient->GetPicture( this->Icon(name), x,y );
}
//______________________________________________________________________________
void NuVldMainFrame::HandleDisabledNeuGEN(void)
{
/*
  // if the neugen interfac is not activated there is no point in:
  if(!fMyConfig->UseNeuGEN()) {
    // ... having the NeuGEN menu enabled
    fMenuNeuGen->DisableEntry(M_NEUGEN_CONFIG_PHYSICS);
    fMenuNeuGen->DisableEntry(M_NEUGEN_CONFIG_PROCESS);
    fMenuNeuGen->DisableEntry(M_NEUGEN_RUN);

    // ... having the fitter tab enabled
    fTabData->SetEnabled(2,kFALSE);
    fMenuFit->DisableEntry(M_FIT_OPEN );
    fMenuFit->DisableEntry(M_FIT_RUN  );
    fMenuFit->DisableEntry(M_FIT_RESET);

    // ... having the neugen buttons enabled
    fNeugenConfigBtn -> SetState(kButtonDisabled);
    fNeugenProcBtn  ->  SetState(kButtonDisabled);
  }
*/
}
//______________________________________________________________________________
void NuVldMainFrame::InitNotify(void)
{
/*
  if(!fMyConfig->UseNeuGEN()) {
      fLog -> AddLine(
             "Dummy GENIE/NeuGEN interface found. NeuGEN options were disabled");
      new GuiMsgBox(gClient->GetRoot(),fMain, 380, 250, kVerticalFrame,
         "Dummy GENIE/NeuGEN interface found. All NeuGEN options were disabled");
  }
*/
}
//______________________________________________________________________________


//______________________________________________________________________________
//--------------------- Window initialization methods --------------------------
//______________________________________________________________________________

//______________________________________________________________________________
void NuVldMainFrame::Init(void)
{
  fDBC  = new DBConnection();
//fNGFP = new NeuGenFitParams();

  fLtxAuth = new TLatex(0.01,0.96, kMajorLabel);
  fLtxLink = new TLatex(0.01,0.92, kMinorLabel);

  fLtxLink -> SetNDC(); // use Normalized Device Coordinates (NDC)
  fLtxLink -> SetTextSize  (0.03);
  fLtxLink -> SetTextColor (9);
  fLtxAuth -> SetNDC(); // use Normalized Device Coordinates (NDC)
  fLtxAuth -> SetTextSize  (0.035);
  fLtxAuth -> SetTextColor (50);
  fLtxAuth -> SetTextFont  (32);

  fPlotterShowIsOn = false;
  fSpline  = 0;

  fNuXSecTab = new GuiNuDataSelectionTab (fMain, fDBC);
  fElXSecTab = new GuiElDataSelectionTab (fDBC);
  fSFTab     = new GuiSFDataSelectionTab (fDBC);
}
//______________________________________________________________________________
void NuVldMainFrame::InitializeHandlers(void)
{
  fHelpHandler    = new GuiHelpHandler    (fMain);
  fDBaseHandler   = new GuiDBHandler      (fMain, fDBC);
  fXmlFileHandler = new GuiXmlFileHandler (fMain, fDBC);

  fStackHandler   = new GuiStackHandler;
//fFitKernel      = new GuiFitKernel;
}
//______________________________________________________________________________
void NuVldMainFrame::InitializeSyslog(void)
{
  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

  syslog->fStatusBar   = fStatusBar;   // init syslog's StatusBar
  syslog->fProgressBar = fProgressBar; //               ProgressBar
  syslog->fLog         = fLog;         //               TextEdit log
}
//______________________________________________________________________________
void NuVldMainFrame::InitializeBrowser(void)
{
  GuiBrowserSingleton * browser = GuiBrowserSingleton::Instance();

  browser->_e_canvas   = fPlotTabEmbCnv;  // init browser's RootEmbeddedCanvas
  browser->_text_edit  = fDataViewer;     //                TextEdit
}
//______________________________________________________________________________
void NuVldMainFrame::ConfigHandlers(void)
{
  fStackHandler -> SetParentMainFrame    (fMain);
  fStackHandler -> AttachDBTableTxtEntry (fStackTableNameTxE);
  fStackHandler -> AttachConfigTxtEntry  (fStackConfigNameTxE);
  fStackHandler -> AttachDBTableCombo    (fTableStackCBx);
  fStackHandler -> AttachConfigCombo     (fConfigStackCBx);
  fStackHandler -> SetDBConnection       (fDBC);
}
//______________________________________________________________________________

//______________________________________________________________________________
//---------------- Decide who is to respond to GUI events ----------------------
//______________________________________________________________________________

//______________________________________________________________________________
void NuVldMainFrame::HandleMenu(Int_t id)
{
  switch(id) {

  case M_FILE_OPEN:             fXmlFileHandler->OpenFile();             break;
  case M_FILE_CLOSE:            fXmlFileHandler->CloseFile();            break;
  case M_FILE_EXIT:             this->CloseWindow();                     break;
  case M_FILE_PARSE:            fXmlFileHandler->ParseFile();            break;
  case M_SAVE_STACKED_DATA:     fStackHandler->SaveStack();              break;
  case M_LOAD_STACKED_DATA:     fStackHandler->LoadStack();              break;
  case M_DATA_CONNECT:          fDBaseHandler->MakeConnection();         break;
  case M_DATA_CLOSE:            fDBaseHandler->CloseConnection();        break;
  case M_DATA_VERIFY:           fDBaseHandler->CheckConnection();        break;
  case M_DATA_DBINFO:           fDBaseHandler->PrintInfo();              break;
  case M_DATA_BOOTSTRAP:        fDBaseHandler->Bootstrap();              break;
  case M_DATA_UPLOAD:           fXmlFileHandler->SendToDBase();          break;
  case M_DATA_QUERY:            fDBaseHandler->QueryWithSqlFromDialog(); break;
  case M_DATA_QUERY_FILE:       fDBaseHandler->QueryWithSqlFromFile();   break;
  case M_DATA_QUERY_DRAW_GUI:   this->DrawDBTable();                     break;
  case M_DATA_QUERY_PRINT_GUI:  this->PrintDBTable();                    break;
  case M_EXPORT_PLOT:           this->HandleSaveCanvas();                break;
  case M_EXPORT_TABLE:          this->HandleSaveTableAsText();           break;
  case M_EXPORT_MODEL:          this->HandleSaveSpline();                break;
  case M_NEUGEN_CONFIG_PHYSICS: this->ConfigNeugenPhysics();             break;
  case M_NEUGEN_CONFIG_PROCESS: this->ConfigNeugenProcess();             break;
  case M_NEUGEN_RUN:            this->RunNeuGen();                       break;
  case M_IMPORT_MODEL:  this->LoadExtXSecPrediction();           break;
  case M_FIT_OPEN:              this->OpenFitterTab();                   break;
  case M_FIT_RUN:               this->RunFitter();                       break;
  case M_FIT_RESET:             this->ResetFitterTab();                  break;
  case M_HELP_ABOUT:            fHelpHandler->NuVldAbout();              break;
  case M_HELP_WWW:              fHelpHandler->NuVldOnline();             break;
  case M_HELP_DURHAM:           fHelpHandler->DurhamOnline();            break;
  case M_HELP_HOWTO_CONN_DBASE: fHelpHandler->HowtoConnDBase();          break;
  case M_HELP_HOWTO_FILL_DBASE: fHelpHandler->HowtoFillDBase();          break;

  case M_VIEW_ENABLE_UNDOCK_MENU:
           fMenuDock->EnableUndock(!fMenuDock->EnableUndock());
           if (fMenuDock->EnableUndock()) {
               fMenuView->CheckEntry(M_VIEW_ENABLE_UNDOCK_MENU);
               fMenuView->EnableEntry(M_VIEW_DOCK_MENU);
               fMenuView->EnableEntry(M_VIEW_UNDOCK_MENU);
           } else {
               fMenuView->UnCheckEntry(M_VIEW_ENABLE_UNDOCK_MENU);
               fMenuView->DisableEntry(M_VIEW_DOCK_MENU);
               fMenuView->DisableEntry(M_VIEW_UNDOCK_MENU);
           }
           break;
  case M_VIEW_ENABLE_UNDOCK_SELECTION_TABS:
           fDBSelDock->EnableUndock(!fDBSelDock->EnableUndock());
           if (fDBSelDock->EnableUndock()) {
               fMenuView->CheckEntry(M_VIEW_ENABLE_UNDOCK_SELECTION_TABS);
               fMenuView->EnableEntry(M_VIEW_DOCK_SELECTION_TABS);
               fMenuView->EnableEntry(M_VIEW_UNDOCK_SELECTION_TABS);
           } else {
               fMenuView->UnCheckEntry(M_VIEW_ENABLE_UNDOCK_SELECTION_TABS);
               fMenuView->DisableEntry(M_VIEW_DOCK_SELECTION_TABS);
               fMenuView->DisableEntry(M_VIEW_UNDOCK_SELECTION_TABS);
           }
           break;
  case M_VIEW_ENABLE_UNDOCK_BOTTOM_FRAME:
           fBottomFrmDock->EnableUndock(!fBottomFrmDock->EnableUndock());
           if (fBottomFrmDock->EnableUndock()) {
               fMenuView->CheckEntry(M_VIEW_ENABLE_UNDOCK_BOTTOM_FRAME);
               fMenuView->EnableEntry(M_VIEW_DOCK_BOTTOM_FRAME);
               fMenuView->EnableEntry(M_VIEW_UNDOCK_BOTTOM_FRAME);
           } else {
               fMenuView->UnCheckEntry(M_VIEW_ENABLE_UNDOCK_BOTTOM_FRAME);
               fMenuView->DisableEntry(M_VIEW_DOCK_BOTTOM_FRAME);
               fMenuView->DisableEntry(M_VIEW_UNDOCK_BOTTOM_FRAME);
           }
           break;

  case M_VIEW_ENABLE_HIDE_MENU:
           fMenuDock->EnableHide(!fMenuDock->EnableHide());
           if (fMenuDock->EnableHide()) {
               fMenuView->CheckEntry(M_VIEW_ENABLE_HIDE_MENU);
               fMenuView->EnableEntry(M_VIEW_HIDE_MENU);
           } else {
               fMenuView->UnCheckEntry(M_VIEW_ENABLE_HIDE_MENU);
               fMenuView->DisableEntry(M_VIEW_HIDE_MENU);
           }
           break;
  case M_VIEW_ENABLE_HIDE_SELECTION_TABS:
           fDBSelDock->EnableHide(!fDBSelDock->EnableHide());
           if (fDBSelDock->EnableHide()) {
               fMenuView->CheckEntry(M_VIEW_ENABLE_HIDE_SELECTION_TABS);
               fMenuView->EnableEntry(M_VIEW_HIDE_SELECTION_TABS);
           } else {
               fMenuView->UnCheckEntry(M_VIEW_ENABLE_HIDE_SELECTION_TABS);
               fMenuView->DisableEntry(M_VIEW_HIDE_SELECTION_TABS);
           }
           break;
  case M_VIEW_ENABLE_HIDE_BOTTOM_FRAME:
           fBottomFrmDock->EnableHide(!fBottomFrmDock->EnableHide());
           if (fBottomFrmDock->EnableHide()) {
               fMenuView->CheckEntry(M_VIEW_ENABLE_HIDE_BOTTOM_FRAME);
               fMenuView->EnableEntry(M_VIEW_HIDE_BOTTOM_FRAME);
           } else {
               fMenuView->UnCheckEntry(M_VIEW_ENABLE_HIDE_BOTTOM_FRAME);
               fMenuView->DisableEntry(M_VIEW_HIDE_BOTTOM_FRAME);
           }
           break;

  case M_VIEW_DOCK_MENU:
          fMenuDock->DockContainer();
          fMenuView->EnableEntry(M_VIEW_UNDOCK_MENU);
          fMenuView->DisableEntry(M_VIEW_DOCK_MENU);
          break;
  case M_VIEW_UNDOCK_MENU:
          fMenuDock->UndockContainer();
          fMenuView->EnableEntry(M_VIEW_DOCK_MENU);
          fMenuView->DisableEntry(M_VIEW_UNDOCK_MENU);
          break;
  case M_VIEW_HIDE_MENU:
          fMenuDock->HideContainer();
          break;

  case M_VIEW_DOCK_SELECTION_TABS:
          fDBSelDock->DockContainer();
          fMenuView->EnableEntry(M_VIEW_UNDOCK_SELECTION_TABS);
          fMenuView->DisableEntry(M_VIEW_DOCK_SELECTION_TABS);
          this->UpdateFrameSize();
          break;
  case M_VIEW_UNDOCK_SELECTION_TABS:
          fDBSelDock->UndockContainer();
          fMenuView->EnableEntry(M_VIEW_DOCK_SELECTION_TABS);
          fMenuView->DisableEntry(M_VIEW_UNDOCK_SELECTION_TABS);
          this->UpdateFrameSize();
          break;
  case M_VIEW_HIDE_SELECTION_TABS:
          fDBSelDock->HideContainer();
          this->UpdateFrameSize();
          break;
  case M_VIEW_SHOW_SELECTION_TABS:
          fDBSelDock->ShowContainer();
          this->UpdateFrameSize();
          break;

  case M_VIEW_DOCK_BOTTOM_FRAME:
          fBottomFrmDock->DockContainer();
          fMenuView->EnableEntry(M_VIEW_UNDOCK_BOTTOM_FRAME);
          fMenuView->DisableEntry(M_VIEW_DOCK_BOTTOM_FRAME);
          this->UpdateFrameSize();
          break;
  case M_VIEW_UNDOCK_BOTTOM_FRAME:
          fBottomFrmDock->UndockContainer();
          fMenuView->EnableEntry(M_VIEW_DOCK_BOTTOM_FRAME);
          fMenuView->DisableEntry(M_VIEW_UNDOCK_BOTTOM_FRAME);
          this->UpdateFrameSize();
          break;
  case M_VIEW_HIDE_BOTTOM_FRAME:
          fBottomFrmDock->HideContainer();
          this->UpdateFrameSize();
          break;
  case M_VIEW_SHOW_BOTTOM_FRAME:
          fBottomFrmDock->ShowContainer();
          this->UpdateFrameSize();
          break;

  default:
       fLog->AddLine( "GUI Event could not be handled" );
       fStatusBar->SetText( "GUI Event could not be handled", 0 );
  }
}
//______________________________________________________________________________
void NuVldMainFrame::UpdateFrameSize(void)
{
  UInt_t dW=0, dH=0;

  bool no_btmfrm = fBottomFrmDock->IsUndocked() || fBottomFrmDock->IsHidden();
  bool no_seltab = fDBSelDock->IsUndocked()     || fDBSelDock->IsHidden();

  if (no_btmfrm && no_seltab) {
     dW = fMainLeftFrame->GetWidth();
     dH = fMainDBottomFrame->GetHeight();
     fMain->Resize(fW0-dW,fH0-dH);
  } else if (no_btmfrm && !no_seltab) {
     dW = fMainLeftFrame->GetWidth();
     dH = fMainDBottomFrame->GetHeight();
     fMain->Resize(fW0,fH0-dH);
  } else {
     fMain->Resize(fW0,fH0);
  }
  gClient->ForceRedraw();
}
//______________________________________________________________________________

//______________________________________________________________________________
//---- Methods for setting / reseting / reading list/combo - box selections ----
//----                   Methods for switching GUI tabs                     ----
//______________________________________________________________________________

//______________________________________________________________________________
void NuVldMainFrame::ResetSqlSelections(void)
{
  // check which SQL tab is active when the reset button is pressed

  if      (fTabSql->GetCurrent() == 0) fNuXSecTab -> ResetSelections();
  else if (fTabSql->GetCurrent() == 1) fElXSecTab -> ResetSelections();
  else if (fTabSql->GetCurrent() == 2) fSFTab     -> ResetSelections();

  this->ResetCommonSelections();
}
//______________________________________________________________________________
void NuVldMainFrame::ClearViewer(void)
{
  if (fTabData->GetCurrent() == 0)
  {
    fPlotTabEmbCnv->GetCanvas()->Clear();
    fPlotTabEmbCnv->GetCanvas()->Update();
    fPlotterShowIsOn = false;
  }
  else if (fTabData->GetCurrent() == 1)
  {
    fDataViewer->Clear();
    fPlotterShowIsOn = false;
  }
}
//______________________________________________________________________________
void NuVldMainFrame::ResetCommonSelections(void)
{
  fShowColorCodeChkB -> SetOn (kTRUE );
  fShowExtLegendChkB -> SetOn (kFALSE);
}
//______________________________________________________________________________
bool NuVldMainFrame::ScaleWithEnergy(void)
{
  string selections = fNuXSecTab->BundleSelectionsInString();

  if( selections.find("scale-with-energy=yes") != string::npos ) return true;
  else return false;
}
//______________________________________________________________________________
string NuVldMainFrame::ErrorOption(void)
{
  string selections = fNuXSecTab->BundleSelectionsInString();

  if      ( selections.find("err-opt=none")      != string::npos ) return "none-noE";
  else if ( selections.find("err-opt=stat+syst") != string::npos ) return "all-noE";
  else if ( selections.find("err-opt=stat")      != string::npos ) return "stat-noE";
  else if ( selections.find("err-opt=syst")      != string::npos ) return "syst-noE";
  else return "";
}
//______________________________________________________________________________
string NuVldMainFrame::PlotVariable(void)
{
  string selections = "";

  if (fTabSql->GetCurrent() == 1)
                      selections = fElXSecTab->BundleSelectionsInString();
  if (fTabSql->GetCurrent() == 2)
                      selections = fSFTab->BundleSelectionsInString();

  vector<string> elements = utils::str::Split(selections,  "$");

  vector<string>::iterator element_iter;

  for(element_iter = elements.begin();
                               element_iter != elements.end(); ++element_iter) {

    if( element_iter->find("DRAW_OPT") != string::npos) {

         vector<string> opt = utils::str::Split(*element_iter,  ";");

         vector<string>::iterator opt_iter;

         for(opt_iter = opt.begin(); opt_iter != opt.end(); ++opt_iter) {

             if(opt_iter->find("plot-var") != string::npos) {

                   vector<string> parts = utils::str::Split(*opt_iter, "=");

                   if(parts.size() == 2) return parts[1];
             }
         }
    }
  }
  return "";
}
//______________________________________________________________________________
void NuVldMainFrame::OpenPlotterTab(void)
{
  fTabData->SetTab(0);
}
//______________________________________________________________________________
void NuVldMainFrame::OpenDataViewerTab(void)
{
  fTabData->SetTab(1);
}
//______________________________________________________________________________
void NuVldMainFrame::OpenFitterTab(void)
{
  fTabData->SetTab(2);
}
//______________________________________________________________________________
void NuVldMainFrame::OpenSessionLogTab(void)
{
  fTabData->SetTab(3);
}
//______________________________________________________________________________

//______________________________________________________________________________
//--------- Methods relevant to the Fitter and the other GUI tabs --------------
//______________________________________________________________________________

//______________________________________________________________________________
void NuVldMainFrame::HandleSaveCanvas(void)
{
  static TString dir(".");

  TGFileInfo fi;
  fi.fFileTypes = kSavedPlotExtensions;
  fi.fIniDir    = StrDup( dir.Data() );

  new TGFileDialog(gClient->GetRoot(), fMain, kFDSave, &fi);

  if( fi.fFilename ) {

     string save_as_filename = string( fi.fFilename );

     ostringstream cmd;
     cmd << "Saving canvas as : " << save_as_filename.c_str();

     fLog->AddLine( cmd.str().c_str() );

     fStatusBar->SetText( cmd.str().c_str(), 0 );
     fStatusBar->SetText( "canvas saved",    1 );

     fPlotTabEmbCnv -> GetCanvas() -> SaveAs( save_as_filename.c_str() );
  }
}
//______________________________________________________________________________
void NuVldMainFrame::HandleSaveSpline(void)
{
  if(!fSpline) {
    fLog       -> AddLine ("No model to export");
    fStatusBar -> SetText ("No model to export", 0 );
    new GuiMsgBox(gClient->GetRoot(),
       fMain, 380, 250, kVerticalFrame, "There is model prediction to export!");
    return;
  }

  static TString dir(".");

  TGFileInfo fi;
  fi.fFileTypes = kSavedSplineExtensions;
  fi.fIniDir    = StrDup( dir.Data() );

  new TGFileDialog(gClient->GetRoot(), fMain, kFDSave, &fi);

  if( fi.fFilename ) {

     string filename = string( fi.fFilename );

     ostringstream cmd;
     cmd << "Saving spline in file : " << filename.c_str();
     fLog       -> AddLine( cmd.str().c_str()    );
     fStatusBar -> SetText( cmd.str().c_str(), 0 );

     if (filename.find(".xml") != string::npos) {
        string splname;
        new GuiTextEntryDialog(gClient->GetRoot(), fMain, 1, 1,
                       kHorizontalFrame, "Type in the spline name", splname);
        if(splname.length() == 0) splname = "NuVldSpline";

        cmd << " with name = " << splname;
        fStatusBar -> SetText( cmd.str().c_str(), 0 );

        fSpline->SaveAsXml(filename,"x","y",splname);
     }

     else if (filename.find(".txt") != string::npos) {
        fSpline->SaveAsText(filename);
     }

     else if (filename.find(".root") != string::npos) {

        string splname;
        new GuiTextEntryDialog(gClient->GetRoot(), fMain, 1, 1,
                          kHorizontalFrame, "Type a spline name", splname);

        cmd << " with name = " << splname;
        fStatusBar -> SetText( cmd.str().c_str(), 0 );

        fSpline->SaveAsROOT(filename, splname, false);
     }
     else {
        new GuiMsgBox(gClient->GetRoot(),
                        fMain, 380, 250, kVerticalFrame, "Unknown file format");
        return;
     }
  }
}
//______________________________________________________________________________
void NuVldMainFrame::HandleSaveTableAsText(void)
{
  GuiBrowserSingleton * browser = GuiBrowserSingleton::Instance();

  this->OpenDataViewerTab();
  browser->TextBrowser()->SaveFile(0,false);
}
//______________________________________________________________________________
void NuVldMainFrame::ResetFitterTab(void)
{
  fXMinNmE -> SetNumber(0);
  fXMaxNmE -> SetNumber(0);

  fFitterCBx->Select(0);

  fFitTxtResults -> Clear();

  fFitTabFuncEmbCnv  -> GetCanvas() -> Clear();
  fFitTabChisqEmbCnv -> GetCanvas() -> Clear();
}
//______________________________________________________________________________

//______________________________________________________________________________
//----- Methods that pop-up windows for data-selection & generator config. -----
//______________________________________________________________________________

//______________________________________________________________________________
void NuVldMainFrame::ConfigNeugenPhysics(void)
{
/*
  fLog       -> AddLine( "Configuring NeuGEN Physics Input"    );
  fStatusBar -> SetText( "Configuring NeuGEN Physics Input", 0 );

  new NeuGenConfigDialog(gClient->GetRoot(), fMain, 400, 200);
*/
}
//______________________________________________________________________________
void NuVldMainFrame::ConfigNeugenProcess(void)
{
/*
  fLog       -> AddLine( "Configuring Modeled Process"    );
  fStatusBar -> SetText( "Configuring Modeled Process", 0 );

  new NeuGenInputDialog(gClient->GetRoot(), fMain, 400, 200);
*/
}
//______________________________________________________________________________
void NuVldMainFrame::RetrieveNeuGenCards(void)
{
/*
  fLog       -> AddLine( "Retrieving NeuGEN Cards"    );
  fStatusBar -> SetText( "Retrieving NeuGEN Cards", 0 );

  NeuGenCards * cards = NeuGenCards::Instance();
  NuVldUserData * user_data = NuVldUserData::Instance();

  if(fUseStackedChkB->GetState() == kButtonDown) {

      if(fConfigStackCBx->GetSelectedEntry()) {

        int entry_id = fConfigStackCBx->GetSelectedEntry()->EntryId();

        if(entry_id >= 0 && entry_id < (int) user_data->NStackedConfigs()) {

           string name = fStackHandler->StackedConfigName(entry_id);

           NGCardPair * cp = user_data -> NeuGenCardPairStack() -> GetCardPair(name);

           cards -> SetConfig ( cp->GetConfig() );
           cards -> SetInputs ( cp->GetInputs() );
        }
      }
  }

  LOG("NuVld", pINFO) << *(cards->CurrInputs());
*/
}
//______________________________________________________________________________
void NuVldMainFrame::RunNeuGen(void)
{
/*
  fLog       -> AddLine( "Running NeuGEN"    );
  fStatusBar -> SetText( "Running NeuGEN", 0 );

  this->RetrieveNeuGenCards();

  NeuGenCards * cards = NeuGenCards::Instance();

  if(this->CheckNeugenCards()) {

    if(fSpline) delete fSpline;

    LOG("NuVld", pDEBUG) << "Starting NeuGEN with configuration : ";
    LOG("NuVld", pDEBUG) << *(cards->CurrConfig());

    // figure out the type of plot
    NGPlotType_t plot_type = cards->CurrInputs()->PlotType();

    // cross section plot
    if(plot_type == e_XSec) {
      int           nbins = cards->CurrInputs()->NBins();
      float         emin  = cards->CurrInputs()->EnergyMin();
      float         emax  = cards->CurrInputs()->EnergyMax();
      NGInteraction intr  = cards->CurrInputs()->GetInteraction();
      NGFinalState  fs    = cards->CurrInputs()->GetFinalState();
      NeuGenCuts    cuts  = cards->CurrInputs()->GetCuts();

      bool is_inclusive = cards->CurrInputs()->Inclusive();

      NeuGenWrapper neugen( cards->CurrConfig() );

      if(is_inclusive)
         fSpline = neugen.XSecSpline( emin, emax, nbins, &intr, &cuts);
      else
         fSpline = neugen.ExclusiveXSecSpline(
                                          emin, emax, nbins, &intr, &fs, &cuts);

      LOG("NuVld", pDEBUG) << "Drawing xsec vs energy ";

      this->DrawSpline (fSpline, fPlotTabEmbCnv);
    }

    // structure function plot
    if(plot_type == e_SF) {

      int           raw_dis_code = cards->CurrInputs()->SFRawDisCode();
      int           nbins        = cards->CurrInputs()->NBins();
      int           A            = cards->CurrInputs()->A();
      float         varmin       = cards->CurrInputs()->PlotVarMin();
      float         varmax       = cards->CurrInputs()->PlotVarMax();
      float         fixvar       = cards->CurrInputs()->SFFixedVar();
      NGSF_t        sftype       = cards->CurrInputs()->SF();
      NGKineVar_t   var          = cards->CurrInputs()->PlotVar();
      NGInteraction intr         = cards->CurrInputs()->GetInteraction();
      NGInitState_t init         = intr.GetInitState();
      NGCcNc_t      ccnc         = intr.GetCCNC();

      LOG("NuVld", pINFO) << "raw dis code = " << raw_dis_code ;

      //NeuGenWrapper neugen( cards->CurrConfig() ); why does this messes with SF's?!!
      NeuGenWrapper neugen;

      fSpline = neugen.StrucFuncSpline(var, varmin, varmax, nbins,
                                 fixvar, A, init, ccnc, sftype, raw_dis_code);

      LOG("NuVld", pDEBUG) << "Drawing structure function vs Q2 or x";

      if(fSpline) this->DrawSpline (fSpline, fPlotTabEmbCnv);
    }

  } else {
      new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                                        "Your NeuGEN cards must be messed up!");
  }
*/
}
//______________________________________________________________________________
void NuVldMainFrame::LoadExtXSecPrediction(void)
{
  fLog       -> AddLine( "Loading xsec prediction from file"    );
  fStatusBar -> SetText( "Loading xsec prediction from file", 0 );

  static TString dir(".");

  const char * kFileExt[] = {"All files", "*", "Ascii files", "*.txt", 0, 0};

  TGFileInfo fi;
  fi.fFileTypes = kFileExt;
  fi.fIniDir    = StrDup(dir.Data());

  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);

  if( fi.fFilename ) {

     string xsec_data_file = string( fi.fFilename );

     ostringstream cmd;
     cmd << "Reading xsec data from file: " << xsec_data_file.c_str();

     fLog       -> AddLine( cmd.str().c_str()   );
     fStatusBar -> SetText( cmd.str().c_str(), 0 );
     fStatusBar -> SetText( "XSec Data File Opened",   1 );

     if(fSpline) delete fSpline;
     fSpline = new Spline(xsec_data_file);

     LOG("NuVld", pDEBUG) << "Drawing xsec vs energy ";

     this->DrawSpline (fSpline, fPlotTabEmbCnv, false);
  }
}
//______________________________________________________________________________
void NuVldMainFrame::DrawSpline(
                   Spline * xs, TRootEmbeddedCanvas * ecanvas, bool show_titles)
{
  bool scale_E = this->ScaleWithEnergy();

  LOG("NuVld", pDEBUG) << "Getting xsec = f (E)";

  TGraph * graph = xs->GetAsTGraph(1000, scale_E);

  graph->SetLineWidth(2);
  graph->SetLineStyle(1);
  graph->SetLineColor(1);

  LOG("NuVld", pDEBUG) << "Checking whether a frame is already drawn";

  if( !fPlotterShowIsOn ) {
      LOG("NuVld", pDEBUG) << "No frame found - Drawing the x,y axes";

      double xmin, xmax, ymin, ymax;

      GraphUtils::range(graph, xmin, xmax, ymin, ymax);

      LOG("NuVld", pDEBUG)
                  << " x = [" << xmin << ", " << xmax << "],"
                  << " y = [" << ymin << ", " << ymax << "] ";

      ecanvas->GetCanvas()->cd();
      ecanvas->GetCanvas()->Clear();

      ecanvas->GetCanvas()->Divide(2,1);

      ecanvas->GetCanvas()->GetPad(1)->SetPad("mplots_pad","",0.01,0.01,0.98,0.99);
      ecanvas->GetCanvas()->GetPad(2)->SetPad("legend_pad","",0.98,0.01,0.99,0.99);

      ecanvas->GetCanvas()->GetPad(1)->SetFillColor(0);
      ecanvas->GetCanvas()->GetPad(1)->SetBorderMode(0);
      ecanvas->GetCanvas()->GetPad(2)->SetFillColor(0);
      ecanvas->GetCanvas()->GetPad(2)->SetBorderMode(0);

      ecanvas->GetCanvas()->GetPad(1)->cd();

      ecanvas->GetCanvas()->GetPad(1)->SetBorderMode(0);

      TH1F * hframe = ecanvas->GetCanvas()->DrawFrame(xmin, TMath::Max(0.,ymin), xmax, 1.2*ymax);

      if(show_titles) {
         hframe->GetXaxis()->SetTitle("Ev (GeV)");
         if(scale_E)
            hframe->GetYaxis()->SetTitle("#sigma/Ev (10^{-38} cm^{2}/GeV");
         else
            hframe->GetYaxis()->SetTitle("#sigma (10^{-38} cm^{2}");
      }
      hframe->Draw();
      graph->Draw("LP");

      fLtxAuth->Draw();
      fLtxLink->Draw();

      fPlotterShowIsOn = true;

      if( xmin > 0 && xmax/xmin > 10. ) ecanvas->GetCanvas()->GetPad(1)->SetLogx();
      if( ymin > 0 && ymax/ymin > 10. ) ecanvas->GetCanvas()->GetPad(1)->SetLogy();
      if( xmin == 0 && xmax > 40.     ) ecanvas->GetCanvas()->GetPad(1)->SetLogx();
      if( ymin == 0 && ymax > 40.     ) ecanvas->GetCanvas()->GetPad(1)->SetLogy();

      ecanvas->GetCanvas()->cd();

  } else {

      LOG("NuVld", pDEBUG) << "Found mplots_pad TPad in TRootEmbeddedCanvas";

      ecanvas->GetCanvas()->GetPad(1)->cd();

      graph->Draw("LP");
      ecanvas->GetCanvas()->GetPad(1)->Update();
/*
      LOG("NuVld", pDEBUG) << "Asking gROOT to finding mplots_pad TPad";

      TPad * pad = (TPad *) gROOT->FindObject("mplots_pad");

      if(pad) {

        LOG("NuVld", pDEBUG)
                  << "Plotting NeuGEN xsecs & updating embedded canvas";

        pad->cd();
        graph->Draw("LP");
        pad->Update();

      } else {
        LOG("NuVld", pERROR) << "gROOT returned a NULL mplots_pad TPad";

        ecanvas->GetCanvas()->GetPad(1)->cd();

        graph->Draw("LP");
        ecanvas->GetCanvas()->GetPad(1)->Update();
      }
*/
  }
  ecanvas->GetCanvas()->Update();
}
//______________________________________________________________________________
bool NuVldMainFrame::CheckNeugenCards(void)
{
/*
// Make sure that what we ask from NeuGEN makes sense

  NeuGenCards * cards = NeuGenCards::Instance();

  bool valid = false;

  int           nbins        = cards->CurrInputs()->NBins();
  int           A            = cards->CurrInputs()->A();
  float         Emin         = cards->CurrInputs()->EnergyMin();
  float         Emax         = cards->CurrInputs()->EnergyMax();
  float         varmin       = cards->CurrInputs()->PlotVarMin();
  float         varmax       = cards->CurrInputs()->PlotVarMax();
  float         fixvar       = cards->CurrInputs()->SFFixedVar();
  NGKineVar_t   var          = cards->CurrInputs()->PlotVar();

  // figure out the type of plot
  NGPlotType_t plot_type = cards->CurrInputs()->PlotType();

  // cross section plot
  if(plot_type == e_XSec) {
    valid = (nbins>2) && (Emin<Emax);
  }
  // structure function plot
  if(plot_type == e_SF) {
    valid = (nbins>2) && (varmin<varmax) && (A>0);
    valid = valid && ((var == e_x   && fixvar > 0) ||
                      (var == e_qqs && fixvar >= 0 && fixvar <= 1));
  }

  return valid;
*/
  return false;
}
//______________________________________________________________________________
DBTable<DBElDiffXSecTableRow> * NuVldMainFrame::FillElDiffXSecTable(void)
{
  fProgressBar->SetPosition(0);

  DBTable<DBElDiffXSecTableRow> * table = 0;

  // read inputs from SQL GUI widgets

  string selections = fElXSecTab->BundleSelectionsInString();

  DBQueryString query_string(selections);

  // create a DBTable loader

  fProgressBar->SetPosition(20);

  if( fDBaseHandler->IsConnected() ) {

     fStatusBar -> SetText("Got connection to SQL Server - Filling DBTable<T>", 1);
     fLog       -> AddLine("Got connection to SQL Server"    );

     DBI dbi( fDBC->SqlServer() );

     // create an empty DBTable

     fProgressBar->SetPosition(40);

     table = new DBTable<DBElDiffXSecTableRow>;

     // fill the table

     dbi.FillTable(table, query_string);

     fProgressBar->SetPosition(80);

  } else {

     fLog -> AddLine("Could not get connection to SQL Server");

     new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250,
                        kVerticalFrame, " No active connection to SQL Server ");
  }
  fProgressBar->SetPosition(100);
  fProgressBar->SetPosition(0);

  return table;
}
//______________________________________________________________________________
DBTable<DBNuXSecTableRow> * NuVldMainFrame::FillNuXSecTable(void)
{
  DBTable<DBNuXSecTableRow> * table = 0;

  // read inputs from SQL GUI widgets

  fProgressBar->SetPosition(10);
  fProgressBar->SetPosition(20);

  if( fDBaseHandler->IsConnected() ) {

    fStatusBar -> SetText( "Found connection to SQL Server", 1 );
    fLog       -> AddLine( "Found connection to SQL Server"    );

    string selections = fNuXSecTab->BundleSelectionsInString();
    fProgressBar->SetPosition(60);

    LOG("NuVld", pDEBUG) << "Selections: " << selections;

    DBQueryString query_string(selections);

    table = new DBTable<DBNuXSecTableRow>;

    DBI dbi( fDBC->SqlServer() );

    dbi.FillTable(table, query_string);

    fProgressBar->SetPosition(80);

  } else {

     fStatusBar -> SetText( "No active connection to SQL Server", 1 );
     fLog       -> AddLine( "No active connection to SQL Server"    );

     fProgressBar->SetPosition(0);

     new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250,
                kVerticalFrame, " No active connection to SQL Server ");
     return 0;
  }
  fProgressBar->SetPosition(100);
  fProgressBar->SetPosition(0);

  //-- done! return the table

  return table;
}
//______________________________________________________________________________
DBTable<DBSFTableRow> * NuVldMainFrame::FillSFTable(void)
{
  DBTable<DBSFTableRow> * table = 0;

  // read inputs from SQL GUI widgets

  fProgressBar->SetPosition(10);
  fProgressBar->SetPosition(20);

  if( fDBaseHandler->IsConnected() ) {

    fStatusBar -> SetText( "Found connection to SQL Server", 1 );
    fLog       -> AddLine( "Found connection to SQL Server"    );

    string selections = fSFTab->BundleSelectionsInString();
    LOG("NuVld", pDEBUG) << "Selections: " << selections;

    fProgressBar->SetPosition(60);

    DBQueryString query_string(selections);

    table = new DBTable<DBSFTableRow>;

    DBI dbi( fDBC->SqlServer() );

    dbi.FillTable(table, query_string);

    fProgressBar->SetPosition(80);

  } else {

     fStatusBar -> SetText( "No active connection to SQL Server", 1 );
     fLog       -> AddLine( "No active connection to SQL Server"    );

     fProgressBar->SetPosition(0);

     new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250,
                         kVerticalFrame, " No active connection to SQL Server ");
     return 0;
  }
  fProgressBar->SetPosition(100);
  fProgressBar->SetPosition(0);

  //-- done! return the table

  return table;
}
//______________________________________________________________________________
void NuVldMainFrame::RetrieveStackedDBTable(void)
{
  LOG("NuVld", pDEBUG) << "Retrieving stacked DBTable<T>";

  NuVldUserData * user_data = NuVldUserData::Instance();

  if(fTableStackCBx->GetSelectedEntry()) {

    int id = fTableStackCBx->GetSelectedEntry()->EntryId();

    LOG("NuVld", pDEBUG) << "Stacked table combo-box id: " << id;

    if(id >= 0 && id < (int) user_data->NStackedDBTables()) {

       string name = fStackHandler->StackedDBTableName(id);

       LOG("NuVld", pDEBUG) << "Stacked table name: " << name;

       // 'status bar' & 'session log' entries

       fStatusBar -> SetText( "Selected stacked data-set:", 0);
       fStatusBar -> SetText(  name.c_str(),                1);

       string log_entry = "Selecting stacked data-set: " + name;

       fLog -> AddLine( log_entry.c_str()  );

       LOG("NuVld", pDEBUG) << "Setting selected stacked table as current";

       user_data->SetStackedDBTableAsCurr(name);

       if( user_data->CurrDBTableIsNull() )
           new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250,
                  kVerticalFrame, " Error! Stacked DBTable was not retrieved. ");
       else {

         DBTableType_t dbtype = user_data->CurrDBTableType();

         LOG("NuVld", pERROR) << "Table type: " << DBTableType::AsString(dbtype);

         if      (dbtype == eDbt_NuXSec    ) fTabSql->SetTab(0);
         else if (dbtype == eDbt_ElDiffXSec) fTabSql->SetTab(1);
         else {
           LOG("NuVld", pERROR) << "Unknown current DBTable<T> type";
         }
       }

    }// entry-id

  } else
      new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                           " You need to select one of the stacked data-sets ");
}
//______________________________________________________________________________
void NuVldMainFrame::SetCurrDBTable(void)
{
   if(fUseStackedChkB->GetState() == kButtonDown) this->RetrieveStackedDBTable();

   else {

      NuVldUserData * user_data = NuVldUserData::Instance();

      if (fTabSql->GetCurrent() == 0)
                           user_data->SetCurrDBTable( this->FillNuXSecTable() );
      if (fTabSql->GetCurrent() == 1)
                       user_data->SetCurrDBTable( this->FillElDiffXSecTable() );
      if (fTabSql->GetCurrent() == 2)
                               user_data->SetCurrDBTable( this->FillSFTable() );
   }
}
//______________________________________________________________________________
void NuVldMainFrame::DrawDBTable(void)
{
   this->OpenPlotterTab();
   this->SetCurrDBTable();

   NuVldUserData * user_data = NuVldUserData::Instance();

   if( !user_data->CurrDBTableIsNull() ) this->DrawCurrentDBTable();
}
//______________________________________________________________________________
void NuVldMainFrame::PrintDBTable(void)
{
   this->OpenDataViewerTab();
   this->SetCurrDBTable();

   NuVldUserData * user_data = NuVldUserData::Instance();

   if( !user_data->CurrDBTableIsNull() ) this->PrintCurrentDBTable();
}
//______________________________________________________________________________
void NuVldMainFrame::DrawCurrentDBTable(void)
{
  fPlotTabEmbCnv->GetCanvas()->cd();

  NuVldUserData * user_data = NuVldUserData::Instance();

  //--- Read selections from the active neutrino SQL GUI tab

  if (fTabSql->GetCurrent() == 0) {

      if( !user_data->CurrDBTableIsNull() ) {

         GuiTableRenderer renderer(fPlotTabEmbCnv);

         renderer.SetScaleWithEnergy( this->ScaleWithEnergy() );
         renderer.SetMultigraph( fShowColorCodeChkB->GetState() == kButtonDown );
         //renderer.SetErrorOption(fNuXSecTab->ReadXSecErrorListbox());
         renderer.SetErrorOption(this->ErrorOption());

         if(fShowExtLegendChkB->GetState() == kButtonDown)
                                       renderer.SetExternalLegend(new TLegend());

         renderer.DrawXSecTable( user_data->NuXSec() );
         renderer.PrintDrawingOptions();

         fLtxAuth->Draw();
         fLtxLink->Draw();

         fPlotterShowIsOn = true;

      } else {

        fStatusBar -> SetText( "pointer to DBTable<T> is null", 1 );
        fLog       -> AddLine( "pointer to DBTable<T> is null"    );

        new GuiMsgBox(gClient->GetRoot(), fMain,
             380, 250, kVerticalFrame, " The table you want to draw is empty ");
      }

  } else

  //--- Read selections from the active electron SQL GUI tab

  if (fTabSql->GetCurrent() == 1) {

      if( !user_data->CurrDBTableIsNull() ) {

         GuiTableRenderer renderer(fPlotTabEmbCnv);

         renderer.SetMultigraph( fShowColorCodeChkB->GetState() == kButtonDown );
//         renderer.SetDrawOption(this->ReadXSecSelectionListbox()); // UPDATE
         renderer.SetPlotVariable(this->PlotVariable());

         if(fShowExtLegendChkB->GetState() == kButtonDown)
                                       renderer.SetExternalLegend(new TLegend());

         renderer.DrawXSecTable( user_data->ElDiffXSec() );
         renderer.PrintDrawingOptions();

         fLtxAuth->Draw();
         fLtxLink->Draw();

         fPlotterShowIsOn = true;

      } else {
        fStatusBar -> SetText( "pointer to DBTable<T> is null", 1 );
        fLog       -> AddLine( "pointer to DBTable<T> is null"    );

        new GuiMsgBox(gClient->GetRoot(), fMain,
             380, 250, kVerticalFrame, " The table you want to draw is empty ");
      }
   } else

  //--- Read selections from the SF SQL GUI tab

  if (fTabSql->GetCurrent() == 2) {

      if( !user_data->CurrDBTableIsNull() ) {

         GuiTableRenderer renderer(fPlotTabEmbCnv);

         renderer.SetMultigraph( fShowColorCodeChkB->GetState() == kButtonDown );
         //renderer.SetDrawOption(this->ReadXSecSelectionListbox());
         renderer.SetPlotVariable(this->PlotVariable());

         if(fShowExtLegendChkB->GetState() == kButtonDown)
                                       renderer.SetExternalLegend(new TLegend());

         renderer.DrawXSecTable( user_data->SF() );
         renderer.PrintDrawingOptions();

         fLtxAuth->Draw();
         fLtxLink->Draw();

         fPlotterShowIsOn = true;

      } else {
        fStatusBar -> SetText( "pointer to DBTable<T> is null", 1 );
        fLog       -> AddLine( "pointer to DBTable<T> is null"    );

        new GuiMsgBox(gClient->GetRoot(), fMain,
             380, 250, kVerticalFrame, " The table you want to draw is empty ");
      }
   }

  fPlotTabEmbCnv->GetCanvas()->Update();
}
//______________________________________________________________________________
void NuVldMainFrame::PrintCurrentDBTable(void)
{
  NuVldUserData * user_data = NuVldUserData::Instance();

  GuiTablePrinter printer;
  printer.ScaleXSecWithEnergy( this->ScaleWithEnergy() );

  if( ! user_data->CurrDBTableIsNull() ) {

     /* neutrino scattering data */
     if (fTabSql->GetCurrent() == 0)
                          printer.PrintTable( user_data->NuXSec() );
     /* electron scattering data */
     else if (fTabSql->GetCurrent() == 1)
                      printer.PrintTable( user_data->ElDiffXSec() );
     /* S/F data */
     else if (fTabSql->GetCurrent() == 2)
                              printer.PrintTable( user_data->SF() );
   } else {
        fStatusBar -> SetText( "pointer to DBTable<T> is null", 1 );
        fLog       -> AddLine( "pointer to DBTable<T> is null"    );
   }
}
//______________________________________________________________________________
void NuVldMainFrame::RunFitter(void)
{
  LOG("NuVld", pDEBUG) << "Running fitter";
/*
  // create a fit-kernel and set GUI fit options

  fFitKernel->Reset();

  fFitKernel->SetFitParams(fNGFP);
  fFitKernel->SetFitRange(fXMinNmE->GetNumber(), fXMaxNmE->GetNumber());
  fFitKernel->SetScaleWithEnergy( this->ScaleWithEnergy() );

  fFitKernel->PrintConfig();

  NuVldUserData * user_data = NuVldUserData::Instance();

  // check whether a fitter was selected

  int    entry_id = fFitterCBx->GetSelectedEntry()->EntryId();
  string fitter   = kFitters[entry_id];

  bool fitter_selected =
       (fFitterCBx->GetSelectedEntry() && (strcmp(fitter.c_str(),"NONE") != 0));

  if(!fitter_selected) {
    new GuiMsgBox(gClient->GetRoot(), fMain,
                  380, 250, kVerticalFrame, " First, you must select a fitter");
    return;
  }

  // if the fit is on stacked data & retrieve them and set as current
  if(fUseStackedChkB->GetState() == kButtonDown) {
    if(fTableStackCBx->GetSelectedEntry()) {

        int entry_id = fTableStackCBx->GetSelectedEntry()->EntryId();

        string xsec_table_name = fStackHandler->StackedDBTableName(entry_id);
        user_data->SetStackedDBTableAsCurr(xsec_table_name);

     }
     else new GuiMsgBox(gClient->GetRoot(), fMain,  380, 250, kVerticalFrame,
                   " If you want to use a stacked data-set, then select one!");
  }


  // have fit data?
  if( user_data->CurrDBTableIsNull() ) {
    new GuiMsgBox(gClient->GetRoot(), fMain,
                 380, 250, kVerticalFrame, " You have selected no data to fit");
    return;
  }

  if (strcmp(fitter.c_str(),"XSEC-NORM") == 0)  fFitKernel->DoSimpleFit(true);
  else {

    int nfitparams = fNGFP->NFittedParams();
    if(nfitparams>0) {
        if      ( strcmp(fitter.c_str(),"SIMPLE") == 0 )     fFitKernel->DoSimpleFit(false);
        else if ( strcmp(fitter.c_str(),"NORM-FLOAT") == 0 ) fFitKernel->DoFloatingNormFit();
        else {
          // should never be here
        }
    }
    else new GuiMsgBox(gClient->GetRoot(), fMain,
           380, 250, kVerticalFrame, " You have selected no parameters to fit");
  }

  this->RunPostFitProcessor();
*/
}
//______________________________________________________________________________
void NuVldMainFrame::RunMcScanner(void)
{
/*
  LOG("NuVld", pDEBUG) << "Running 2-D scanner";

  // create a fit-kernel ans set GUI fit options

  fFitKernel->Reset();

  fFitKernel->SetFitParams(fNGFP);
  fFitKernel->SetFitRange(fXMinNmE->GetNumber(), fXMaxNmE->GetNumber());
  fFitKernel->SetScaleWithEnergy( this->ScaleWithEnergy() );

  fFitKernel->PrintConfig();

  // run the scanner & plot

  bool ok = fFitKernel->MCParamScanning();

  if(ok) {
    LOG("NuVld", pDEBUG) << "Plotting xsec boundaries";

    this->PlotXSecBoundaries( fPlotTabEmbCnv->GetCanvas(), false );

    fPlotTabEmbCnv->GetCanvas()->GetPad(1)->Update();
    fPlotTabEmbCnv->GetCanvas()->Update();

    this->OpenPlotterTab();

  } else new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
             " You have not configured the multi-D MC param scanner properly ");
*/
}
//______________________________________________________________________________
void NuVldMainFrame::Run2dScanner(void)
{
/*
  LOG("NuVld", pDEBUG) << "Running MC scanner";

  // create a fit-kernel ans set GUI fit options

  fFitKernel->Reset();

  fFitKernel->SetFitParams(fNGFP);
  fFitKernel->SetFitRange(fXMinNmE->GetNumber(), fXMaxNmE->GetNumber());
  fFitKernel->SetScaleWithEnergy( this->ScaleWithEnergy() );

  fFitKernel->PrintConfig();

  // run the scanner & plot

  bool ok = fFitKernel->ChisqScan2D();

  if(ok) {
    fFitTabFuncEmbCnv->GetCanvas()->cd();
    fFitTabFuncEmbCnv->GetCanvas()->Clear();

    gStyle->SetPalette(1);
    fFitKernel->chisq2d->Draw("COLZ");

    fFitTabFuncEmbCnv->GetCanvas()->Update();

  } else new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                           " You have not configured the 2d-scanner properly ");
*/
}
//______________________________________________________________________________
void NuVldMainFrame::Run1dScanner(void)
{
/*
  LOG("NuVld", pDEBUG) << "Running 1-D scanner";

  // create a fit-kernel ans set GUI fit options

  fFitKernel->Reset();

  fFitKernel->SetFitParams(fNGFP);
  fFitKernel->SetFitRange(fXMinNmE->GetNumber(), fXMaxNmE->GetNumber());
  fFitKernel->SetScaleWithEnergy( this->ScaleWithEnergy() );

  fFitKernel->PrintConfig();

  // run the scanner & plot

  fFitKernel->MCParamScanning();

  LOG("NuVld", pDEBUG) << "Plotting xsec boundaries";

  this->PlotXSecBoundaries( fFitTabFuncEmbCnv->GetCanvas(), true );

  fFitTabFuncEmbCnv->GetCanvas()->Update();

  fFitTabChisqEmbCnv->GetCanvas()->cd();
  fFitTabChisqEmbCnv->GetCanvas()->Clear();

  fFitKernel->ChisqScan1D();
  //double xmin = fFitKernel->chisq1d->GetX()[TMath::LocMin(fFitKernel->chisq1d->GetN(), fFitKernel->chisq1d->GetX())];
  //double xmax = fFitKernel->chisq1d->GetX()[TMath::LocMax(fFitKernel->chisq1d->GetN(), fFitKernel->chisq1d->GetX())];
  //double ymin = fFitKernel->chisq1d->GetY()[TMath::LocMin(fFitKernel->chisq1d->GetN(), fFitKernel->chisq1d->GetY())];
  //double ymax = fFitKernel->chisq1d->GetY()[TMath::LocMax(fFitKernel->chisq1d->GetN(), fFitKernel->chisq1d->GetY())];

  //TH1F * hframe = fFitTabFuncEmbCnv->GetCanvas()->DrawFrame(xmin,ymin,xmax,ymax);

  //hframe->Draw();

  fFitKernel->chisq1d->Draw("ALP");

  fLtxAuth->Draw();
  fLtxLink->Draw();

  fFitTabChisqEmbCnv->GetCanvas()->Update();
*/
}
//______________________________________________________________________________
void NuVldMainFrame::PlotXSecBoundaries(TCanvas * /*c*/, bool /*clear*/)
{
/*
  fFitKernel->lowb->SetLineWidth(2);
  fFitKernel->highb->SetLineWidth(2);
  fFitKernel->lowb->SetLineStyle(2);
  fFitKernel->highb->SetLineStyle(2);

  bool new_pad = false;

  if(! c->GetPad(1) || clear) {

    LOG("NuVld", pWARN) << "NULL pad";

    new_pad = true;

    c->Clear();
    c->Divide(2,1);

    c->GetPad(1)->SetPad("mplots_pad","",0.01,0.01,0.98,0.99);
    c->GetPad(2)->SetPad("legend_pad","",0.98,0.01,0.99,0.99);

    c->GetPad(1)->SetFillColor(0);
    c->GetPad(1)->SetBorderMode(0);
    c->GetPad(2)->SetFillColor(0);
    c->GetPad(2)->SetBorderMode(0);
  }
  c->GetPad(1)->cd();

  if( new_pad ) {

      LOG("NuVld", pDEBUG) << "No frame found - Drawing the x,y axes";

      gStyle->SetOptTitle(0);

      double xmin = fFitKernel->lowb->GetX()[TMath::LocMin(fFitKernel->lowb->GetN(), fFitKernel->lowb->GetX())];
      double xmax = fFitKernel->highb->GetX()[TMath::LocMax(fFitKernel->highb->GetN(), fFitKernel->highb->GetX())];
      double ymin = fFitKernel->lowb->GetY()[TMath::LocMin(fFitKernel->lowb->GetN(), fFitKernel->lowb->GetY())];
      double ymax = fFitKernel->highb->GetY()[TMath::LocMax(fFitKernel->highb->GetN(), fFitKernel->highb->GetY())];

      TH1F * hframe = c->GetPad(1)->DrawFrame(xmin,TMath::Max(0.,ymin),xmax,ymax);

      hframe->Draw();

      fFitKernel->lowb->Draw("LP");
      fFitKernel->highb->Draw("LP");

      if( xmin > 0 && xmax/xmin > 10. ) c->GetPad(1)->SetLogx();
      if( ymin > 0 && ymax/ymin > 10. ) c->GetPad(1)->SetLogy();
      if( xmin == 0 && xmax > 40. ) c->GetPad(1)->SetLogx();
      if( ymin == 0 && ymax > 40. ) c->GetPad(1)->SetLogy();

      fLtxAuth->Draw();
      fLtxLink->Draw();

      c->GetPad(1)->Update();

  } else {

      fFitKernel->lowb->Draw("LP");
      fFitKernel->highb->Draw("LP");
  }
*/
}
//______________________________________________________________________________
void NuVldMainFrame::RunPostFitProcessor(void)
{
/*
// After the fitting is being done, this methods updates all the GUI elements
// waiting to show the fitting results

  LOG("NuVld", pDEBUG) << "Running post fit processor";

  NuVldUserData * user_data = NuVldUserData::Instance();

  LOG("NuVld", pDEBUG)<< "Drawing fitted DBTable<T> in fPlotTabEmbCnv";

//  fPlotTabEmbCnv->GetCanvas()->cd();
//  fPlotTabEmbCnv->GetCanvas()->GetPad(1)->cd();
  this->DrawDBTable();

  LOG("NuVld", pDEBUG) << "Getting & drawing function in fPlotTabEmbCnv";

  fPlotTabEmbCnv->GetCanvas()->cd(1);

  TF1 * func = fFitKernel -> FitFunction();

  if (func) {

    func->Draw("same");

    fPlotTabEmbCnv->GetCanvas()->cd(1)->Update();
    fPlotTabEmbCnv->GetCanvas()->Update();

    LOG("NuVld", pDEBUG)
            << "Drawing fitted DBTable<T> & fit function in fFitTabFuncEmbCnv";

    fFitTabFuncEmbCnv->GetCanvas()->cd();

    GuiTableRenderer renderer(fFitTabFuncEmbCnv);

    renderer.SetScaleWithEnergy( this->ScaleWithEnergy() );
    renderer.SetMultigraph(false);
    //renderer.SetErrorOption(fNuXSecTab->ReadXSecErrorListbox());
    renderer.SetErrorOption(this->ErrorOption());
    renderer.DrawXSecTable( user_data->NuXSec() );
    renderer.PrintDrawingOptions();

    fFitTabFuncEmbCnv->GetCanvas()->cd(1);

    func->Draw("same");

    fLtxAuth->Draw();

    fFitTabFuncEmbCnv->GetCanvas()->GetPad(1)->Update();
    fFitTabFuncEmbCnv->GetCanvas()->Update();

  } else {
    LOG("NuVld", pERROR) << "Fit function is NULL";
  }

  this->PrintFitParameters();
  this->DrawResiduals();
*/
}
//______________________________________________________________________________
void NuVldMainFrame::PrintFitParameters(void)
{
/*
  LOG("NuVld", pDEBUG) << "Getting fit func & and printing fit params";

  TF1 * func = fFitKernel -> FitFunction();
  if (func) {
    fFitTxtResults->AddLine("-----------------------------------------------");
    fFitTxtResults->AddLine("PHYSICS PARAMS:");
    for(int i = 0; i < kNNGFitParams; i++) {
       ostringstream fitparam;

       if( fNGFP->IsFitted(i) ) {
          fitparam << fNGFP->ParamAsString(i) << " = "
                   << func->GetParameter(i) << " +/- " << func->GetParError(i);
       } else {
          fitparam << fNGFP->ParamAsString(i) << " = "
                   << func->GetParameter(i) << " **FIXED**";
       }
       fFitTxtResults->AddLine(fitparam.str().c_str());
    }
    fFitTxtResults->AddLine("XSEC NORM:");
    ostringstream norm;
    norm << func->GetParameter(kNNGFitParams)
                      << " +/- " << func->GetParError(kNNGFitParams);
    fFitTxtResults->AddLine(norm.str().c_str());

    fFitTxtResults->AddLine(Concat("* x^2     = ", func->GetChisquare()));
    fFitTxtResults->AddLine(Concat("* ndf     = ", func->GetNDF()));

    if(func->GetNDF() > 0) {
       float chisq_ndf = func->GetChisquare() / func->GetNDF();
       fFitTxtResults->AddLine(Concat("* x^2/ndf = ", chisq_ndf));
    }

    fFitTxtResults->AddLine(Concat("* Prob    = ", func->GetProb()));

    fFitTxtResults->AddLine("-----------------------------------------------");

  } else {
    LOG("NuVld", pERROR) << "Fit function is NULL";
  }
*/
}
//______________________________________________________________________________
void NuVldMainFrame::DrawResiduals(void)
{
/*
  // get residuals as graph

  LOG("NuVld", pDEBUG) << "Getting residuals as graph";

  TGraph * gr = fFitKernel->GetResidualsAsGraph();

  if(gr) {
    // draw graph frame

    LOG("NuVld", pDEBUG) << "Drawing graph frame";

    double xmin =  999999, ymin =  999999;
    double xmax = -999999, ymax = -999999;

    for(int i = 0; i < gr->GetN() ; i++) {

      xmin = TMath::Min(xmin, gr->GetX()[i]);
      xmax = TMath::Max(xmax, gr->GetX()[i]);
      ymin = TMath::Min(ymin, gr->GetY()[i]);
      ymax = TMath::Max(ymax, gr->GetY()[i]);
    }
    fFitTabChisqEmbCnv->GetCanvas()->cd();

    TH1F * hframe = (TH1F *)
            fFitTabChisqEmbCnv->GetCanvas()->DrawFrame(xmin, ymin, xmax, ymax);

    hframe->Draw();

    if(xmin>0 && xmax/xmin > 10) fFitTabChisqEmbCnv->GetCanvas()->SetLogx();

    LOG("NuVld", pDEBUG) << "Drawing graph & updating the embedded canvas";

    gr->SetMarkerStyle(3);
    gr->SetMarkerSize(1);
    gr->Draw("P");

    fFitTabChisqEmbCnv->GetCanvas()->Update();

    delete gr;
  }
*/
}
//______________________________________________________________________________
