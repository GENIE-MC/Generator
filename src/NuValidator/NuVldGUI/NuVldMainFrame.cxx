//_____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldMainFrame

\brief    NuValidator GUI prototype - main frame

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <TStyle.h>
#include <TFile.h>
#include <TSQLResult.h>
#include <TLatex.h>
#include <TH2F.h>

#include "DBUtils/DBI.h"
#include "DBUtils/DBQueryString.h"
#include "DBUtils/SqlUtils.hh"
#include "Messenger/Messenger.h"
#include "Facades/NeuGenWrapper.h"
#include "NuVldGUI/NuVldMainFrame.h"
#include "NuVldGUI/NuVldUserData.h"
#include "NuVldGUI/GuiTablePrinter.h"
#include "NuVldGUI/GuiTableRenderer.h"
#include "NuVldGUI/MsgBox.h"
#include "NuVldGUI/MultiLineMsgBox.h"
#include "NuVldGUI/NeuGenFitParamsDialog.h"
#include "NuVldGUI/YNQuestionBox.h"
#include "NuVldGUI/TextEntryDialog.h"
#include "NuVldGUI/NeuGenConfigDialog.h"
#include "NuVldGUI/NeuGenInputDialog.h"
#include "NuVldGUI/NeuGenCards.h"
#include "NuVldGUI/DBConnectionDialog.h"
#include "NuVldGUI/SysLogSingleton.h"
#include "NuVldGUI/BrowserSingleton.h"
#include "NuVldGUI/GraphUtils.hh"
#include "NuVldGUI/vDataSelectionDialog.h"
#include "NuVldGUI/vMeasurementListDialog.h"
#include "NuVldGUI/NuVldConstants.h"
#include "Utils/StringUtils.h"
#include "Utils/GUIUtils.h"

#ifndef _NUVLD_DB_UTILS_H_
#include "NuVldGUI/DBUtils.hh"
#endif

using std::string;
using std::ostringstream;
using std::vector;

using namespace genie;
using namespace genie::string_utils;
using namespace genie::nuvld;
using namespace genie::nuvld::constants;

ClassImp(NuVldMainFrame)

//______________________________________________________________________________
NuVldMainFrame::NuVldMainFrame(const TGWindow * p, UInt_t w, UInt_t h) :
TGMainFrame(p, w, h)
{
  fMain = new TGMainFrame(p,w,h);
  fMain->Connect("CloseWindow()",
                         "genie::nuvld::NuVldMainFrame", this, "CloseWindow()");

  this->Init();
  this->InitializeHandlers();    // initialize GUI event handlers

  //-- define TGLayoutHints

  this->DefineLayoutHints();

  //-- add menu bar

  fMenu = this->BuildMenuBar();

  fMain->AddFrame(fMenu, fMenuBarLt);

  //-- instantiate main frames (below menu & above the status bar)

  fMainFrame       = new TGCompositeFrame(fMain,       1, 1, kVerticalFrame);
  fMainTopFrame    = new TGCompositeFrame(fMainFrame, 3, 3, kHorizontalFrame);
  fMainBottomFrame = new TGCompositeFrame(fMainFrame, 3, 3, kHorizontalFrame);

  //-- TOP FRAME: add image buttons frame

  fImgBtnGrpFrm = this->BuildUpperButtonFrame();

  fMainTopFrame -> AddFrame( fImgBtnGrpFrm  );

  //-- BOOTOM FRAME:
  //         left  : neutrino / electron scattering data SQL tabs
  //         right : plotter / data viewer / session log tabs

  fMainLeftFrame   = new TGCompositeFrame(fMainBottomFrame, 3, 3, kVerticalFrame);
  fMainRightFrame  = new TGCompositeFrame(fMainBottomFrame, 3, 3, kVerticalFrame);

  fMainBottomFrame  -> AddFrame ( fMainLeftFrame,    fMLeftFrameLt  );
  fMainBottomFrame  -> AddFrame ( fMainRightFrame,   fMRightFrameLt );

  fTabSql  = this->BuildSqlTab();
  fTabData = this->BuildDataTab();

  fMainRightFrame -> AddFrame ( fTabData, fDataTabLt );
  fMainLeftFrame  -> AddFrame ( fTabSql,  fSqlTabLt  );

  this->AddCommonCheckButtons();

  //-- add Selection-Stacking GUI widgets

  fStackHFrm = this->BuildSelectionStackFrame();

  fMainRightFrame->AddFrame( fStackHFrm, fSelStackLt );

  //-- add Progress Bar & lower row of image buttons

  fProgressBarHFrm = this->BuildLowerButtonFrame();

  fMainRightFrame->AddFrame( fProgressBarHFrm, fProgressBarLt );

  //-- add top/bottom main frames to main frame, and main frame to main

  fMainFrame -> AddFrame ( fMainTopFrame    );
  fMainFrame -> AddFrame ( fMainBottomFrame );

  fMain -> AddFrame ( fMainFrame  );

  //-- add Status Bar

  fStatusBar = this->BuildStatusBar();

  fMain->AddFrame(fStatusBar, fStatusBarLt);

  //-- initialize

  this->ResetSqlSelections();    // set all selection widgets to default values
  this->ResetFitterTab();        // set all selections to default values
  
  this->InitializeSyslog();      // initialize "system-log" singleton
  this->InitializeBrowser();     // initialize "data-browser" singleton
  this->ConfigHandlers();  

  fMain->SetWindowName("GENIE/NuValidator");

  fMain->MapSubwindows();

  fMain->Resize( fMain->GetDefaultSize() );

  fMain->MapWindow();
}
//______________________________________________________________________________
NuVldMainFrame::~NuVldMainFrame()
{
  fMain->Cleanup();
  delete fMain;
}
//______________________________________________________________________________

//______________________________________________________________________________
//----------------- Methods for GUI window construction ------------------------
//______________________________________________________________________________

//______________________________________________________________________________
void NuVldMainFrame::DefineLayoutHints(void)
{
  ULong_t hintMenuBarLayout       = kLHintsTop | kLHintsLeft | kLHintsExpandX;
  ULong_t hintMenuBarItemLayout   = kLHintsTop | kLHintsLeft;
  ULong_t hintMenuBarHelpLayout   = kLHintsTop | kLHintsRight;
  ULong_t hintTabPlotterLayout    = kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintTabFitterLayout     = kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintTabLogLayout        = kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintTabDataViewerLayout = kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintTabNuSqlLayout      = kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintTabElSqlLayout      = kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintTabDataLayout       = kLHintsTop | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintTabSqlLayout        = kLHintsTop | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintExitButtonLayout    = kLHintsRight;
  ULong_t hintLeftButtonsLayout   = kLHintsLeft;
  ULong_t hintStatusBarLayout     = kLHintsBottom | kLHintsLeft | kLHintsExpandX;
  ULong_t hintMLeftFrameLayout    = kLHintsCenterY;
  ULong_t hintMRightFrameLayout   = kLHintsTop | kLHintsExpandX | kLHintsExpandY;
  ULong_t hintProgressBarLayout   = kLHintsCenterX;
  ULong_t hintSelectStackLayout   = kLHintsCenterX;
  ULong_t hintFitLeftFrameLayout  = kLHintsTop | kLHintsLeft;
  ULong_t hintFitRightFrameLayout = kLHintsTop | kLHintsRight;

  fMenuBarLt       = new TGLayoutHints(hintMenuBarLayout,       0, 0,  1, 1);
  fMenuBarItemLt   = new TGLayoutHints(hintMenuBarItemLayout,   0, 4,  0, 0);
  fMenuBarHelpLt   = new TGLayoutHints(hintMenuBarHelpLayout);
  fPlotterTabLt    = new TGLayoutHints(hintTabPlotterLayout,    5, 5, 10, 1);
  fLogTabLt        = new TGLayoutHints(hintTabLogLayout,        5, 5, 10, 1);
  fFitterTabLt     = new TGLayoutHints(hintTabFitterLayout,     5, 5, 10, 1);
  fDataViewTabLt   = new TGLayoutHints(hintTabDataViewerLayout, 5, 5, 10, 1);
  fNuSqlTabLt      = new TGLayoutHints(hintTabNuSqlLayout,      5, 5, 10, 1);
  fElSqlTabLt      = new TGLayoutHints(hintTabElSqlLayout,      5, 5, 10, 1);
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
  TGMenuBar * menu_bar = new TGMenuBar(fMain, 1, 1, kHorizontalFrame);

  // Menu: File

  fMenuFile = new TGPopupMenu(gClient->GetRoot());

  fMenuFile->AddEntry ("&Open XML",               M_FILE_OPEN  );
  fMenuFile->AddEntry ("&Close XML",              M_FILE_CLOSE );
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry ("&Test XML File Parsing",  M_FILE_PARSE );
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry ("&Load Stack",             M_LOAD_STACKED_DATA );
  fMenuFile->AddEntry ("&Save Stack",             M_SAVE_STACKED_DATA );
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry ("E&xit",                   M_FILE_EXIT  );

  fMenuFile->DisableEntry(M_FILE_CLOSE); // disable "close" until next "open"

  fMenuFile->Connect("Activated(Int_t)",
                                  "genie::nuvld::NuVldMainFrame",
                                                     this,"HandleMenu(Int_t)");

  fXmlFileHandler->AttachFileMenu(fMenuFile); // to (de/)activate open/close

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

  // Menu: NeuGEN

  fMenuNeuGen = new TGPopupMenu(gClient->GetRoot());

  fMenuNeuGen->AddEntry("Physics",       M_NEUGEN_CONFIG_PHYSICS);
  fMenuNeuGen->AddEntry("Process",       M_NEUGEN_CONFIG_PROCESS);
  fMenuNeuGen->AddEntry("Run",           M_NEUGEN_RUN);
  fMenuNeuGen->AddSeparator();
  fMenuNeuGen->AddEntry("Load external", M_NEUGEN_LOAD_EXTERNAL);

  fMenuNeuGen->Connect("Activated(Int_t)",
                                   "genie::nuvld::NuVldMainFrame",
                                                    this,"HandleMenu(Int_t)");
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

  menu_bar = new TGMenuBar(fMain, 1, 1, kHorizontalFrame);

  menu_bar->AddPopup("&File",     fMenuFile,   fMenuBarItemLt);
  menu_bar->AddPopup("&Database", fMenuDBase,  fMenuBarItemLt);
  menu_bar->AddPopup("&NeuGEN",   fMenuNeuGen, fMenuBarItemLt);
  menu_bar->AddPopup("&Fit",      fMenuFit,    fMenuBarItemLt);
  menu_bar->AddPopup("&Help",     fMenuHelp,   fMenuBarHelpLt);

  return menu_bar;
}
//______________________________________________________________________________
TGGroupFrame * NuVldMainFrame::BuildUpperButtonFrame(void)
{
  TGGroupFrame * grpf = new TGGroupFrame(fMainTopFrame, "", kVerticalFrame);

  fBtnMatrixLt = new TGMatrixLayout(grpf, 0, 15, 1);

  grpf->SetLayoutManager( fBtnMatrixLt );

  this->CreateUpperFrameButtons(grpf); // create buttons
  this->SetUpperFrameButtonText();            // set tool-tip text
  this->ConnectUpperFrameButtons();           // connect Click() with a method

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

  TGTab * tab = new TGTab(fMainLeftFrame, 1, 1);

  //-- tab: SQL GUI widgets for v scattering data

  tf = tab->AddTab( "vN" );

  fTabNuSql = new TGCompositeFrame(tf, width, height, kVerticalFrame);

  this->FillNuSqlFrame();

  tf->AddFrame( fTabNuSql, fNuSqlTabLt );

  //-- tab: SQL GUI widgets for e scattering data

  tf = tab->AddTab( "eN" );

  fTabElSql = new TGCompositeFrame(tf, width, height, kVerticalFrame);

  this->FillElSqlFrame();

  tf -> AddFrame( fTabElSql, fElSqlTabLt );

  return tab;
}
//______________________________________________________________________________
void NuVldMainFrame::FillNuSqlFrame(void)
{
  fNuXSecErrGrpFrm   = new TGGroupFrame(fTabNuSql, "Cross Section Err",  kVerticalFrame);
  fNuExpGrpFrm       = new TGGroupFrame(fTabNuSql, "Experiment",         kVerticalFrame);
  fNuXSecGrpFrm      = new TGGroupFrame(fTabNuSql, "Cross Section",      kVerticalFrame);
  fEnergyGrpFrm      = new TGGroupFrame(fTabNuSql, "Energy Range (GeV)", kVerticalFrame);
  fNuInitStateGrpFrm = new TGGroupFrame(fTabNuSql, "Initial State",      kVerticalFrame);

  fEnergyMatrixLt = new TGMatrixLayout(fEnergyGrpFrm, 0, 2, 2);
  fEnergyGrpFrm->SetLayoutManager( fEnergyMatrixLt );

  fNuXSecErrLBx = new TGListBox(fNuXSecErrGrpFrm,   2);
  fNuExpLBx     = new TGListBox(fNuExpGrpFrm,       2);
  fNuProcLBx    = new TGListBox(fNuXSecGrpFrm,      2);
  fNuTypeLBx    = new TGListBox(fNuInitStateGrpFrm, 2);
  fNuTgtLBx     = new TGListBox(fNuInitStateGrpFrm, 2);

  gui_utils::FillListBox( fNuXSecErrLBx, kXSecErrType    );
  gui_utils::FillListBox( fNuExpLBx,     kExperimentName );
  gui_utils::FillListBox( fNuProcLBx,    kProcName       );
  gui_utils::FillListBox( fNuTypeLBx,    kNuType         );
  gui_utils::FillListBox( fNuTgtLBx,     kTarget         );

  fNuXSecErrLBx -> Resize (100,  60);
  fNuExpLBx     -> Resize (100,  60);
  fNuProcLBx    -> Resize (100,  60);
  fNuTypeLBx    -> Resize (100,  60);
  fNuTgtLBx     -> Resize (100,  60);

  fNuXSecErrLBx -> SetMultipleSelections( false );
  fNuExpLBx     -> SetMultipleSelections( true  );
  fNuProcLBx    -> SetMultipleSelections( true  );
  fNuTypeLBx    -> SetMultipleSelections( true  );
  fNuTgtLBx     -> SetMultipleSelections( true  );

  fAllNuExpChkB    = new TGCheckButton(fNuExpGrpFrm,       "Select all", 71);
  fAllNuProcChkB   = new TGCheckButton(fNuXSecGrpFrm,      "Select all", 72);
  fAllNuTypesChkB  = new TGCheckButton(fNuInitStateGrpFrm, "Select all", 73);
  fAllNuTgtChkB    = new TGCheckButton(fNuInitStateGrpFrm, "Select all", 74);

  fAllNuExpChkB    -> Connect("Clicked()","genie::nuvld::NuVldMainFrame",
                                                       this,"SelectAllNuExp()");
  fAllNuProcChkB   -> Connect("Clicked()","genie::nuvld::NuVldMainFrame",
                                                      this,"SelectAllNuXSec()");
  fAllNuTypesChkB  -> Connect("Clicked()","genie::nuvld::NuVldMainFrame",
                                                    this,"SelectAllNuProbes()");
  fAllNuTgtChkB    -> Connect("Clicked()","genie::nuvld::NuVldMainFrame",
                                                   this,"SelectAllNuTargets()");

  fNuXSecErrGrpFrm   -> AddFrame( fNuXSecErrLBx   );
  fNuExpGrpFrm       -> AddFrame( fNuExpLBx       );
  fNuExpGrpFrm       -> AddFrame( fAllNuExpChkB   );
  fNuXSecGrpFrm      -> AddFrame( fNuProcLBx      );
  fNuXSecGrpFrm      -> AddFrame( fAllNuProcChkB  );
  fNuInitStateGrpFrm -> AddFrame( fNuTypeLBx      );
  fNuInitStateGrpFrm -> AddFrame( fAllNuTypesChkB );
  fNuInitStateGrpFrm -> AddFrame( fNuTgtLBx       );
  fNuInitStateGrpFrm -> AddFrame( fAllNuTgtChkB   );

  fScaleWithEvChkB  = new TGCheckButton(fTabNuSql, "Scale With Energy", 75);

  fEMinNmE = new TGNumberEntry(
                   fEnergyGrpFrm, kEmin, 6, 1, TGNumberFormat::kNESReal);
  fEMaxNmE = new TGNumberEntry(
                   fEnergyGrpFrm, kEmax, 6, 1, TGNumberFormat::kNESReal);

  fMinELb = new TGLabel(fEnergyGrpFrm, new TGString( "min:"));
  fMaxELb = new TGLabel(fEnergyGrpFrm, new TGString( "max:"));

  fEnergyGrpFrm -> AddFrame ( fMinELb  );
  fEnergyGrpFrm -> AddFrame ( fEMinNmE );
  fEnergyGrpFrm -> AddFrame ( fMaxELb  );
  fEnergyGrpFrm -> AddFrame ( fEMaxNmE );

  fNuTabBtnSpacerLb = new TGLabel(fTabNuSql, new TGString(" "));

  fShowFullNuDialogTBtn   = new TGTextButton (fTabNuSql, "More data selections... ", 76);
  fShowExpertNuDialogTBtn = new TGTextButton (fTabNuSql, "Select by dbase entry...", 77);

  fShowFullNuDialogTBtn->Connect("Clicked()", "genie::nuvld::NuVldMainFrame",
                                       this, "PopupNuDataSelectionDialog()");

  fShowExpertNuDialogTBtn->Connect("Clicked()", "genie::nuvld::NuVldMainFrame",
                                       this, "PopupNuMeasurementListDialog()");

  //-- bottom/left side: add all parent frames

  fTabNuSql -> AddFrame( fNuXSecErrGrpFrm        );
  fTabNuSql -> AddFrame( fNuExpGrpFrm            );
  fTabNuSql -> AddFrame( fNuXSecGrpFrm           );
  fTabNuSql -> AddFrame( fEnergyGrpFrm           );
  fTabNuSql -> AddFrame( fNuInitStateGrpFrm      );
  fTabNuSql -> AddFrame( fNuTabBtnSpacerLb       );
  fTabNuSql -> AddFrame( fShowFullNuDialogTBtn   );
  fTabNuSql -> AddFrame( fShowExpertNuDialogTBtn );
  fTabNuSql -> AddFrame( fNuTabBtnSpacerLb       );
  fTabNuSql -> AddFrame( fScaleWithEvChkB        );
}
//______________________________________________________________________________
void NuVldMainFrame::FillElSqlFrame(void)
{
  fElExpGrpFrame = new TGGroupFrame(fTabElSql, "Experiment", kVerticalFrame);
  fElTgGrpFrm    = new TGGroupFrame(fTabElSql, "Target",     kVerticalFrame);

  fElExpLBx = new TGListBox(fElExpGrpFrame,  222);
  fElTgtLBx = new TGListBox(fElTgGrpFrm,  223);

  gui_utils::FillListBox( fElExpLBx,  kElExperiment );
  gui_utils::FillListBox( fElTgtLBx,  kElTarget     );

  fElExpLBx -> Resize (100,  60);
  fElTgtLBx -> Resize (100,  50);

  fElExpLBx -> SetMultipleSelections( true );
  fElTgtLBx -> SetMultipleSelections( true );

  fAllElExpChkB  = new TGCheckButton(fElExpGrpFrame, "Select all", 401);
  fAllElTgtChkB  = new TGCheckButton(fElTgGrpFrm, "Select all", 402);

  fAllElExpChkB  -> Connect("Clicked()",
                                 "genie::nuvld::NuVldMainFrame", this,"SelectAllElExp()");
  fAllElTgtChkB  -> Connect("Clicked()",
                             "genie::nuvld::NuVldMainFrame", this,"SelectAllElTargets()");

  fElExpGrpFrame -> AddFrame( fElExpLBx     );
  fElExpGrpFrame -> AddFrame( fAllElExpChkB );

  fElTgGrpFrm -> AddFrame( fElTgtLBx     );
  fElTgGrpFrm -> AddFrame( fAllElTgtChkB );

  fTabElSql -> AddFrame( fElExpGrpFrame );
  fTabElSql -> AddFrame( fElTgGrpFrm    );

  // build all (E,Ep,W2,Q2,Theta,v) min/max selections

  for(int iframe = 0; iframe < kNElVarRangeFrames; iframe++) {

     fElVarRangeGrpFrm.push_back( new TGGroupFrame(
                  fTabElSql, kElVarFrameName[iframe], kVerticalFrame) );

     fElVarRangeLt.push_back( new TGMatrixLayout(
                                  fElVarRangeGrpFrm[iframe], 0, 2, 1) );

     fElVarRangeGrpFrm[iframe]->SetLayoutManager(
                                                 fElVarRangeLt[iframe] );


     fElVarMinNmEV.push_back( new TGNumberEntry(
                  fElVarRangeGrpFrm[iframe], kElVarMin[iframe],
                                             6, 1, TGNumberFormat::kNESReal) );

     fElVarMaxNmEV.push_back( new TGNumberEntry(
                  fElVarRangeGrpFrm[iframe], kElVarMax[iframe],
                                             6, 1, TGNumberFormat::kNESReal) );

     fElVarRangeGrpFrm[iframe] -> AddFrame ( fElVarMinNmEV   [iframe] );
     fElVarRangeGrpFrm[iframe] -> AddFrame ( fElVarMaxNmEV   [iframe] );

     fTabElSql -> AddFrame ( fElVarRangeGrpFrm[iframe] );
  }

  // x-variable

  fElDrawXGrpFrm = new TGGroupFrame(fTabElSql,"x-variable:", kVerticalFrame);

  fElDrawXCBx = new TGComboBox(fElDrawXGrpFrm, 412);

  gui_utils::FillComboBox( fElDrawXCBx, kElVarFrameName );

  fElDrawXCBx -> Resize (115, 20);

  fElDrawXGrpFrm -> AddFrame (fElDrawXCBx);
  fTabElSql      -> AddFrame (fElDrawXGrpFrm);

  // Init state

  fAllElExpChkB -> SetOn (kTRUE);
  fAllElTgtChkB -> SetOn (kTRUE);

  fElDrawXCBx->Select(5);
  
  this->SelectAllElExp();  
  this->SelectAllElTargets();
}
//______________________________________________________________________________
void NuVldMainFrame::AddCommonCheckButtons(void)
{
  fShowColorCodeChkB = new TGCheckButton(fMainLeftFrame, "Color-Code plot",  72);
  fShowExtLegendChkB = new TGCheckButton(fMainLeftFrame, "External legend",  73);
  fUseStackedChkB    = new TGCheckButton(fMainLeftFrame, "Use stacked",      74);

  fMainLeftFrame -> AddFrame( fShowColorCodeChkB );
  fMainLeftFrame -> AddFrame( fShowExtLegendChkB );
  fMainLeftFrame -> AddFrame( fUseStackedChkB    );
}
//______________________________________________________________________________
TGTab * NuVldMainFrame::BuildDataTab(void)
{
// Build tabs

  TGCompositeFrame * tf = 0;

  unsigned int width  = 650;
  unsigned int height = 550;

  TGTab * tab = new TGTab(fMainRightFrame, 1, 1);

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

  gui_utils::FillComboBox( fFitterCBx, kFitters );

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

  fDoFitBtn     -> SetToolTipText( "Do the fit",                               1);
  fPrmScanBtn   -> SetToolTipText( "multi-D parameter space MC scanning",      2);
  fPrmScan1dBtn -> SetToolTipText( "1-D parameter space scanning",             3);
  fPrmScan2dBtn -> SetToolTipText( "2-D parameter space scanning",             4);
  fResetFitBtn  -> SetToolTipText( "Reset fitter selections" ,                 5);

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
  new NeuGenFitParamsDialog(
                    gClient->GetRoot(), fMain, 380, 250, kVerticalFrame, fNGFP);
}
//______________________________________________________________________________
TGHorizontalFrame * NuVldMainFrame::BuildSelectionStackFrame(void)
{
  TGHorizontalFrame * hf = new TGHorizontalFrame(fMainRightFrame, 10, 10);

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
  TGHorizontalFrame * hf = new TGHorizontalFrame(fMainRightFrame, 10, 10);

  // lower frame: left hand-side picture buttons

  fSelResetBtn  = new TGPictureButton(hf, this->Pic("reset",      32,32));
  fViewClearBtn = new TGPictureButton(hf, this->Pic("viewclear",  32,32));
  fDrawDataBtn  = new TGPictureButton(hf, this->Pic("draw_data",  32,32));
  fPrintDataBtn = new TGPictureButton(hf, this->Pic("print_data", 32,32));
  fNeugenRunBtn = new TGPictureButton(hf, this->Pic("neugen_run", 32,32));
  fSaveBtn      = new TGPictureButton(hf, this->Pic("save",       32,32));

  fSelResetBtn  -> SetToolTipText( "Reset selection" , 1);
  fViewClearBtn -> SetToolTipText( "Clear view" , 1);
  fDrawDataBtn  -> SetToolTipText( "Query dbase for selected experiments & draw data", 1);
  fPrintDataBtn -> SetToolTipText( "Query dbase for selected experiments & print data in text format", 1);
  fNeugenRunBtn -> SetToolTipText( "Run NeuGEN with selected inputs & draw its prediction", 1);
  fSaveBtn      -> SetToolTipText( "Save graph", 1);

  fSelResetBtn  -> Connect("Clicked()", "genie::nuvld::NuVldMainFrame",
                                                 this,"ResetSqlSelections()");
  fViewClearBtn -> Connect("Clicked()", "genie::nuvld::NuVldMainFrame",
                                                        this,"ClearViewer()");
  fDrawDataBtn  -> Connect("Clicked()", "genie::nuvld::NuVldMainFrame",
                                                      this,"DrawDBTable()");
  fPrintDataBtn -> Connect("Clicked()", "genie::nuvld::NuVldMainFrame",
                                                     this,"PrintDBTable()");
  fNeugenRunBtn -> Connect("Clicked()", "genie::nuvld::NuVldMainFrame",
                                                           this,"RunNulook()");
  fSaveBtn      -> Connect("Clicked()", "genie::nuvld::NuVldMainFrame",
                                                   this,"HandleSaveCanvas()");

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

//______________________________________________________________________________
//--------------------- Window initialization methods --------------------------
//______________________________________________________________________________

//______________________________________________________________________________
void NuVldMainFrame::Init(void)
{
  fDBC  = new DBConnection();
  fNGFP = new NeuGenFitParams();
  
  fLtxAuth = new TLatex(0.01,0.96,
              "GENIE - C.Andreopoulos (CCLRC,Rutherford), H.Gallagher (Tufts)");
  fLtxLink = new TLatex(0.01,0.92,
                          "http://hepunx.rl.ac.uk/~candreop/generators/GENIE/");
  
  fLtxLink -> SetNDC(); // use Normalized Device Coordinates (NDC)
  fLtxLink -> SetTextSize  (0.03);
  fLtxLink -> SetTextColor (9);
  fLtxAuth -> SetNDC(); // use Normalized Device Coordinates (NDC)
  fLtxAuth -> SetTextSize  (0.035);
  fLtxAuth -> SetTextColor (50);
  fLtxAuth -> SetTextFont  (32);

  fPlotterShowIsOn = false;
  _xsec_vs_energy  = 0;

  _v_slc_dialog_requires_attn = false;
  _active_v_slc_dialog        = 0;
}
//______________________________________________________________________________
void NuVldMainFrame::InitializeHandlers(void)
{
  fHelpHandler    = new GuiHelpHandler    (fMain);
  fDBaseHandler   = new GuiDBHandler      (fMain, fDBC);
  fXmlFileHandler = new GuiXmlFileHandler (fMain, fDBC);

  fStackHandler   = new GuiStackHandler;
  fFitKernel      = new GuiFitKernel;
}
//______________________________________________________________________________
void NuVldMainFrame::InitializeSyslog(void)
{
  SysLogSingleton * syslog = SysLogSingleton::Instance();

  syslog->fStatusBar   = fStatusBar;   // init syslog's StatusBar
  syslog->fProgressBar = fProgressBar; //               ProgressBar
  syslog->fLog         = fLog;         //               TextEdit log
}
//______________________________________________________________________________
void NuVldMainFrame::InitializeBrowser(void)
{
  BrowserSingleton * browser = BrowserSingleton::Instance();

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
  case M_NEUGEN_CONFIG_PHYSICS: this->ConfigNeugenPhysics();             break;
  case M_NEUGEN_CONFIG_PROCESS: this->ConfigNeugenProcess();             break;
  case M_NEUGEN_RUN:            this->RunNulook();                       break;
  case M_NEUGEN_LOAD_EXTERNAL:  this->LoadExtXSecPrediction();           break;
  case M_FIT_OPEN:              this->OpenFitterTab();                   break;
  case M_FIT_RUN:               this->RunFitter();                       break;
  case M_FIT_RESET:             this->ResetFitterTab();                  break;
  case M_HELP_ABOUT:            fHelpHandler->NuVldAbout();              break;
  case M_HELP_WWW:              fHelpHandler->NuVldOnline();             break;
  case M_HELP_DURHAM:           fHelpHandler->DurhamOnline();            break;
  case M_HELP_HOWTO_CONN_DBASE: fHelpHandler->HowtoConnDBase();          break;
  case M_HELP_HOWTO_FILL_DBASE: fHelpHandler->HowtoFillDBase();          break;

   default:
       fLog->AddLine( "GUI Event could not be handled" );
       fStatusBar->SetText( "GUI Event could not be handled", 0 );
  }
}
//______________________________________________________________________________

//______________________________________________________________________________
//---- Methods for setting / reseting / reading list/combo - box selections ----
//----                   Methods for switching GUI tabs                     ----
//______________________________________________________________________________

//______________________________________________________________________________
void NuVldMainFrame::SelectAllNuExp(void)
{
  if(fAllNuExpChkB->GetState() == kButtonDown)
                                  gui_utils::SelectAllListBoxEntries(fNuExpLBx);
  else gui_utils::ResetAllListBoxSelections(fNuExpLBx);

  fNuExpLBx->SelectionChanged();

  gClient->NeedRedraw(fNuExpLBx->GetContainer());
}
//______________________________________________________________________________
void NuVldMainFrame::SelectAllNuXSec(void)
{
  if(fAllNuProcChkB->GetState() == kButtonDown)
                                 gui_utils::SelectAllListBoxEntries(fNuProcLBx);
  else gui_utils::ResetAllListBoxSelections(fNuProcLBx);

  fNuProcLBx->SelectionChanged();

  gClient->NeedRedraw(fNuProcLBx->GetContainer());
}
//______________________________________________________________________________
void NuVldMainFrame::SelectAllNuProbes(void)
{
  if(fAllNuTypesChkB->GetState() == kButtonDown)
                                 gui_utils::SelectAllListBoxEntries(fNuTypeLBx);
  else gui_utils::ResetAllListBoxSelections(fNuTypeLBx);

  fNuTypeLBx->SelectionChanged();

  gClient->NeedRedraw(fNuTypeLBx->GetContainer());
}
//______________________________________________________________________________
void NuVldMainFrame::SelectAllNuTargets(void)
{
  if(fAllNuTgtChkB->GetState() == kButtonDown)
                                  gui_utils::SelectAllListBoxEntries(fNuTgtLBx);
  else gui_utils::ResetAllListBoxSelections(fNuTgtLBx);

  fNuTgtLBx->SelectionChanged();

  gClient->NeedRedraw(fNuTgtLBx->GetContainer());
}
//______________________________________________________________________________
void NuVldMainFrame::SelectAllElExp(void)
{
  if(fAllElExpChkB->GetState() == kButtonDown)
                                  gui_utils::SelectAllListBoxEntries(fElExpLBx);
  else gui_utils::ResetAllListBoxSelections(fElExpLBx);

  fElExpLBx->SelectionChanged();

  gClient->NeedRedraw(fElExpLBx->GetContainer());
}
//______________________________________________________________________________
void NuVldMainFrame::SelectAllElTargets(void)
{
  if(fAllElTgtChkB->GetState() == kButtonDown)
                                 gui_utils::SelectAllListBoxEntries(fElTgtLBx);
  else gui_utils::ResetAllListBoxSelections(fElTgtLBx);

  fElTgtLBx->SelectionChanged();

  gClient->NeedRedraw(fElTgtLBx->GetContainer());
}
//______________________________________________________________________________
void NuVldMainFrame::ResetSqlSelections(void)
{
  // check which SQL tab is active when the reset button is pressed

  if      (fTabSql->GetCurrent() == 0) this->ResetNuSqlSelections();
  else if (fTabSql->GetCurrent() == 1) this->ResetElSqlSelections();

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
void NuVldMainFrame::ResetNuSqlSelections(void)
{
  gui_utils::ResetAllListBoxSelections( fNuExpLBx  );
  gui_utils::ResetAllListBoxSelections( fNuProcLBx );
  gui_utils::ResetAllListBoxSelections( fNuTypeLBx );
  gui_utils::ResetAllListBoxSelections( fNuTgtLBx  );

  fEMinNmE->SetNumber(kEmin);
  fEMaxNmE->SetNumber(kEmax);

  fNuXSecErrLBx  -> Select (2);

  fAllNuExpChkB   -> SetOn (kTRUE);
  fAllNuProcChkB  -> SetOn (kTRUE);
  fAllNuTypesChkB -> SetOn (kTRUE);
  fAllNuTgtChkB   -> SetOn (kTRUE);

  this->SelectAllNuExp();
  this->SelectAllNuXSec();
  this->SelectAllNuProbes();
  this->SelectAllNuTargets();
}
//______________________________________________________________________________
void NuVldMainFrame::ResetElSqlSelections(void)
{
  gui_utils::ResetAllListBoxSelections( fElExpLBx );
  gui_utils::ResetAllListBoxSelections( fElTgtLBx );

  for(int iframe = 0; iframe < kNElVarRangeFrames; iframe++) {

     fElVarMinNmEV[iframe] -> SetNumber ( kElVarMin[iframe] );
     fElVarMaxNmEV[iframe] -> SetNumber ( kElVarMax[iframe] );
  }

  fAllElExpChkB -> SetOn (kTRUE);
  fAllElTgtChkB -> SetOn (kTRUE);

  fElDrawXCBx -> Select (5);

  this->SelectAllElExp();
  this->SelectAllElTargets();
}
//______________________________________________________________________________
void NuVldMainFrame::ResetCommonSelections(void)
{
  fScaleWithEvChkB   -> SetOn (kTRUE );
  fShowColorCodeChkB -> SetOn (kTRUE );
  fShowExtLegendChkB -> SetOn (kFALSE);
}
//______________________________________________________________________________
string NuVldMainFrame::NuDataSelections(void)
{
  string selections = "";

  if(_v_slc_dialog_requires_attn)
       selections = _active_v_slc_dialog->BundleSelectionsInString();
  else
       selections = this->NuTabBundleSelectionsInString();

  return selections;
}
//______________________________________________________________________________
bool NuVldMainFrame::ScaleWithEnergy(void)
{
  NuVldUserData * user_data = NuVldUserData::Instance();

  if( user_data->CurrDBTableIsNull() ) 
                       return (fScaleWithEvChkB->GetState() == kButtonDown);  
  else {
     string selections = this->NuDataSelections();

     if( selections.find("scale-with-energy") != string::npos ) return true;
     else return false;
  }
}
//______________________________________________________________________________
string NuVldMainFrame::PlotVariable(void)
{
  string selections = this->ElDataSelections();

  vector<string> elements = ParserUtils::split(selections,  "$");

  vector<string>::iterator element_iter;

  for(element_iter = elements.begin();
                               element_iter != elements.end(); ++element_iter) {

    if( element_iter->find("DRAW_OPT") != string::npos) {

         vector<string> opt = ParserUtils::split(*element_iter,  ";");

         vector<string>::iterator opt_iter;

         for(opt_iter = opt.begin(); opt_iter != opt.end(); ++opt_iter) {

             if(opt_iter->find("plot-var") != string::npos) {

                   vector<string> parts = ParserUtils::split(*opt_iter, "=");

                   if(parts.size() == 2) return parts[1];
             }
         }
    }
  }
  return "";
}
//______________________________________________________________________________
string NuVldMainFrame::ReadXSecSelectionListbox(void)
{
  ostringstream err;

  TGLBEntry * selected_entry = fNuXSecErrLBx->GetSelectedEntry();

  if(selected_entry) {

    err << kXSecErrDrawOpt[ selected_entry->EntryId() ] << "-noE";

    fLog->AddLine( Concat(
     "Cross Section Errors - List Box selection: ", selected_entry->EntryId()) );

  } else {

    err << "allXsec-noE";
    fLog->AddLine( "No Cross Section Error Selection - setting default" );
  }

  LOG("NuVld", pDEBUG) << "error selection = " << err.str().c_str();

  return err.str();
}
//______________________________________________________________________________
string NuVldMainFrame::ElDataSelections(void)
{
  string selections = this->ElTabBundleSelectionsInString();

  return selections;
}
//______________________________________________________________________________
string NuVldMainFrame::NuTabBundleSelectionsInString(void)
{
  ostringstream options;

  options << "KEY-LIST:" << this->NuTabBundleKeyListInString()  << "$"
          << "CUTS:"     << this->NuTabBundleCutsInString()     << "$"
          << "DRAW_OPT:" << this->NuTabBundleDrawOptInString()  << "$"
          << "DB-TYPE:vN-XSec";

  return options.str();
}
//______________________________________________________________________________
string NuVldMainFrame::NuTabBundleKeyListInString(void)
{
  // Read experiment name selections
  string experiments = gui_utils::ListBoxSelectionAsString(
                                                fNuExpLBx, kExperimentMySQLName);
  // Read xsec selections
  string xsecs = gui_utils::ListBoxSelectionAsString(fNuProcLBx, kProcMySQLName);
  // Read neutrino selections
  string nus = gui_utils::ListBoxSelectionAsString(fNuTypeLBx, kNuTypeMySQLName);
  // Read target selections
  string targets = gui_utils::ListBoxSelectionAsString(
                                                    fNuTgtLBx, kTargetMySQLName);

  fLog->AddLine( Concat("requested experiments : ", experiments.c_str()) );
  fLog->AddLine( Concat("requested measurements : ", xsecs.c_str()) );
  fLog->AddLine( Concat("requested neutrino beams : ", nus.c_str()) );
  fLog->AddLine( Concat("requested targets : ", targets.c_str()) );

  // Build key list
  string key_list = SqlUtils::build_v_key_list(
                fDBC->SqlServer(), experiments, xsecs, nus, targets);

  return key_list;
}
//______________________________________________________________________________
string NuVldMainFrame::NuTabBundleCutsInString(void)
{
  float Emin = fEMinNmE->GetNumber();
  float Emax = fEMaxNmE->GetNumber();

  ostringstream cuts;

  cuts << "Emin=" << Emin << ";" << "Emax=" << Emax;

  return cuts.str();
}
//______________________________________________________________________________
string NuVldMainFrame::NuTabBundleDrawOptInString(void)
{
  if(fScaleWithEvChkB->GetState() == kButtonDown) return "scale-with-energy";
  else return "";
}
//______________________________________________________________________________
string NuVldMainFrame::ElTabBundleSelectionsInString(void)
{
  ostringstream options;

  options << "KEY-LIST:" << this->ElTabBundleKeyListInString() << "$"
          << "CUTS:"     << this->ElTabBundleCutsInString()    << "$"
          << "DRAW_OPT:" << this->ElTabBundleDrawOptInString() << "$"
          << "DB-TYPE:eN-Diff-XSec";
          
  return options.str();
}
//______________________________________________________________________________
string NuVldMainFrame::ElTabBundleKeyListInString(void)
{
  // Read experiment name selections
  string exprm = gui_utils::ListBoxSelectionAsString(fElExpLBx, kElExperiment);

  // Read target selections
  string targets = gui_utils::ListBoxSelectionAsString(fElTgtLBx, kElTarget);

  // Build key list
  string key_list = SqlUtils::build_e_key_list(fDBC->SqlServer(), exprm, targets);

  return key_list;
}
//______________________________________________________________________________
string NuVldMainFrame::ElTabBundleCutsInString(void)
{
  float E_min       =  fElVarMinNmEV[0]->GetNumber();
  float E_max       =  fElVarMaxNmEV[0]->GetNumber();
  float EP_min      =  fElVarMinNmEV[1]->GetNumber();
  float EP_max      =  fElVarMaxNmEV[1]->GetNumber();
  float Theta_min   =  fElVarMinNmEV[2]->GetNumber();
  float Theta_max   =  fElVarMaxNmEV[2]->GetNumber();
  float Q2_min      =  fElVarMinNmEV[3]->GetNumber();
  float Q2_max      =  fElVarMaxNmEV[3]->GetNumber();
  float W2_min      =  fElVarMinNmEV[4]->GetNumber();
  float W2_max      =  fElVarMaxNmEV[4]->GetNumber();
  float Nu_min      =  fElVarMinNmEV[5]->GetNumber();
  float Nu_max      =  fElVarMaxNmEV[5]->GetNumber();
  float Epsilon_min =  fElVarMinNmEV[6]->GetNumber();
  float Epsilon_max =  fElVarMaxNmEV[6]->GetNumber();
  float Gamma_min   =  fElVarMinNmEV[7]->GetNumber();
  float Gamma_max   =  fElVarMaxNmEV[7]->GetNumber();
  float x_min       =  fElVarMinNmEV[8]->GetNumber();
  float x_max       =  fElVarMaxNmEV[8]->GetNumber();

  ostringstream cuts;

  cuts << "E_min="       << E_min       << ";"
       << "E_max="       << E_max       << ";"
       << "EP_min="      << EP_min      << ";"
       << "EP_max="      << EP_max      << ";"
       << "Theta_min="   << Theta_min   << ";"
       << "Theta_max="   << Theta_max   << ";"
       << "Q2_min="      << Q2_min      << ";"
       << "Q2_max="      << Q2_max      << ";"
       << "W2_min="      << W2_min      << ";"
       << "W2_max="      << W2_max      << ";"
       << "Nu_min="      << Nu_min      << ";"
       << "Nu_max="      << Nu_max      << ";"
       << "Epsilon_min=" << Epsilon_min << ";"
       << "Epsilon_max=" << Epsilon_max << ";"
       << "Gamma_min="   << Gamma_min   << ";"
       << "Gamma_max="   << Gamma_max   << ";"
       << "x_min="       << x_min       << ";"
       << "x_max="       << x_max       << ";";

  return cuts.str();
}
//______________________________________________________________________________
string NuVldMainFrame::ElTabBundleDrawOptInString(void)
{
  int selected = fElDrawXCBx->GetSelected();

  const char * plot_var = kElVarMySQLName[selected];

  ostringstream draw_opt;

  draw_opt << "plot-var=" << plot_var;

  return draw_opt.str();
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
  fLog       -> AddLine( "Configuring NeuGEN Physics Input"    );
  fStatusBar -> SetText( "Configuring NeuGEN Physics Input", 0 );

  new NeuGenConfigDialog(gClient->GetRoot(), fMain, 400, 200);
}
//______________________________________________________________________________
void NuVldMainFrame::ConfigNeugenProcess(void)
{
  fLog       -> AddLine( "Configuring Modeled Process"    );
  fStatusBar -> SetText( "Configuring Modeled Process", 0 );

  new NeuGenInputDialog(gClient->GetRoot(), fMain, 400, 200);
}
//______________________________________________________________________________
void NuVldMainFrame::PopupNuDataSelectionDialog(void)
{
  bool IsConnected;

  if( !fDBC->SqlServer() ) IsConnected = false;
  else IsConnected = fDBC->SqlServer()->IsConnected();

  if(IsConnected) {

     if(!_v_slc_dialog_requires_attn) {

      _active_v_slc_dialog = new vDataSelectionDialog(
                  gClient->GetRoot(), fMain, _v_slc_dialog_requires_attn,
                                   750, 500, kHorizontalFrame, fDBC );
     } else {

       new MsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
           "Another selection dialog has locked my attention. Close it first.");
     }

  } else {
      new MsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                   "You must be connected to the data-base to use this option");
  }
}
//______________________________________________________________________________
void NuVldMainFrame::PopupNuMeasurementListDialog(void)
{
  bool IsConnected;

  if( !fDBC->SqlServer() ) IsConnected = false;
  else IsConnected = fDBC->SqlServer()->IsConnected();

  if(IsConnected) {

     if(!_v_slc_dialog_requires_attn) {

      _active_v_slc_dialog = new vMeasurementListDialog(
                  gClient->GetRoot(), fMain, _v_slc_dialog_requires_attn,
                                               650, 400, kVerticalFrame, fDBC );
     } else {

       new MsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
           "Another selection dialog has locked my attention. Close it first.");
     }

  } else {
      new MsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                   "You must be connected to the data-base to use this option");
  }
}
//______________________________________________________________________________
void NuVldMainFrame::RunNulook(void)
{
  fLog       -> AddLine( "Running NeuGEN"    );
  fStatusBar -> SetText( "Running NeuGEN", 0 );

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

  if(this->CheckNeugenCards()) {

    if(_xsec_vs_energy) delete _xsec_vs_energy;

    LOG("NuVld", pDEBUG) << "NeuGEN with configuration : ";
    LOG("NuVld", pDEBUG) << *(cards->CurrConfig());

    NeuGenWrapper neugen( cards->CurrConfig() );

    int           nbins = cards->CurrInputs()->NBins();
    float         emin  = cards->CurrInputs()->EnergyMin();
    float         emax  = cards->CurrInputs()->EnergyMax();
    NGInteraction intr  = cards->CurrInputs()->GetInteraction();
    NGFinalState  fs    = cards->CurrInputs()->GetFinalState();
    NeuGenCuts    cuts  = cards->CurrInputs()->GetCuts();

    bool is_inclusive = cards->CurrInputs()->Inclusive();

    if(is_inclusive)
        _xsec_vs_energy = neugen.XSecSpline( emin, emax, nbins, &intr, &cuts);
    else
        _xsec_vs_energy = neugen.ExclusiveXSecSpline(
                                          emin, emax, nbins, &intr, &fs, &cuts);

    LOG("NuVld", pDEBUG) << "Drawing xsec vs energy ";

    this->DrawNeugenXSecVsEnergy (_xsec_vs_energy, fPlotTabEmbCnv);

  } else {
      new MsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                                        "Your NeuGEN cards must be messed up!");
  }    
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

     if(_xsec_vs_energy) delete _xsec_vs_energy;
     
     _xsec_vs_energy = new XSecVsEnergy;
     _xsec_vs_energy->LoadFromFile( xsec_data_file.c_str() );
     
     LOG("NuVld", pDEBUG) << "Drawing xsec vs energy ";

     this->DrawNeugenXSecVsEnergy (_xsec_vs_energy, fPlotTabEmbCnv, false);
  }
}
//______________________________________________________________________________
void NuVldMainFrame::DrawNeugenXSecVsEnergy(
             XSecVsEnergy * xs, TRootEmbeddedCanvas * ecanvas, bool show_titles)
{
  bool scale_E = this->ScaleWithEnergy();

  LOG("NuVld", pDEBUG) << "Getting xsec = f (E)";
  
  TGraph * graph = xs->GetAsGraph(1000, scale_E);

  graph->SetLineWidth(2);
  graph->SetLineStyle(1);
  graph->SetLineColor(1);

  LOG("NuVld", pDEBUG) << "Checking whether a frame is already drawn";
  
  //NuVldUserData * user_data = NuVldUserData::Instance();

  //if( user_data->CurrDBTableIsNull() ) {
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
// Make sure that what we ask from NeuGEN makes sense

  NeuGenCards * cards = NeuGenCards::Instance();

  bool vld_NBins  = cards->CurrInputs()->NBins() > 0;
  bool vld_Erange = cards->CurrInputs()->EnergyMax() >
                                               cards->CurrInputs()->EnergyMin();
 
  bool valid = vld_NBins && vld_Erange;

  return valid;
}
//______________________________________________________________________________
DBTable<eDiffXSecTableRow> * NuVldMainFrame::FillElDiffXSecTable(void)
{
  fProgressBar->SetPosition(0);

  DBTable<eDiffXSecTableRow> * table = 0;

  // read inputs from SQL GUI widgets

  string selections = ElTabBundleSelectionsInString();

  DBQueryString query_string(selections);
  
  // create a DBTable loader

  fProgressBar->SetPosition(20);

  if( fDBaseHandler->IsConnected() ) {

     fStatusBar -> SetText("Got connection to SQL Server - Filling DBTable<T>", 1);
     fLog       -> AddLine("Got connection to SQL Server"    );

     DBI dbi( fDBC->SqlServer() );

     // create an empty DBTable

     fProgressBar->SetPosition(40);

     table = new DBTable<eDiffXSecTableRow>;

     // fill the table

     dbi.FillTable(table, query_string);

     fProgressBar->SetPosition(80);

  } else {

     fLog -> AddLine("Could not get connection to SQL Server");

     new MsgBox(gClient->GetRoot(), fMain, 380, 250,
                        kVerticalFrame, " No active connection to SQL Server ");
  }

  fProgressBar->SetPosition(100);
  fProgressBar->SetPosition(0);

  return table;
}
//______________________________________________________________________________
DBTable<vXSecTableRow> * NuVldMainFrame::FillNuXSecTable(void)
{
  DBTable<vXSecTableRow> * table = 0;

  // read inputs from SQL GUI widgets

  fProgressBar->SetPosition(10);
  fProgressBar->SetPosition(20);

  if( fDBaseHandler->IsConnected() ) {

    fStatusBar -> SetText( "Found connection to SQL Server", 1 );
    fLog       -> AddLine( "Found connection to SQL Server"    );

    string selections = "";

    if(_v_slc_dialog_requires_attn)
    {
       selections = _active_v_slc_dialog->BundleSelectionsInString();
       fProgressBar->SetPosition(60);
    }
    else
    {
       selections = this->NuTabBundleSelectionsInString();
       fProgressBar->SetPosition(60);
    }

    LOG("NuVld", pDEBUG) << "Selections: " << selections;

    DBQueryString query_string(selections);

    table = new DBTable<vXSecTableRow>;

    DBI dbi( fDBC->SqlServer() );
    
    dbi.FillTable(table, query_string);
    
    fProgressBar->SetPosition(80);

  } else {

     fStatusBar -> SetText( "No active connection to SQL Server", 1 );
     fLog       -> AddLine( "No active connection to SQL Server"    );

     fProgressBar->SetPosition(0);

     new MsgBox(gClient->GetRoot(), fMain, 380, 250,
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
           new MsgBox(gClient->GetRoot(), fMain, 380, 250,
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
      new MsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
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
         renderer.SetErrorOption(this->ReadXSecSelectionListbox());

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

        new MsgBox(gClient->GetRoot(), fMain,
             380, 250, kVerticalFrame, " The table you want to draw is empty ");
      }

  } else

  //--- Read selections from the active electron SQL GUI tab

  if (fTabSql->GetCurrent() == 1) {

      if( !user_data->CurrDBTableIsNull() ) {

         GuiTableRenderer renderer(fPlotTabEmbCnv);

         renderer.SetMultigraph( fShowColorCodeChkB->GetState() == kButtonDown );
         renderer.SetDrawOption(this->ReadXSecSelectionListbox());
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

        new MsgBox(gClient->GetRoot(), fMain,
             380, 250, kVerticalFrame, " The table you want to draw is empty ");
      }
   }//e

  fPlotTabEmbCnv->GetCanvas()->Update();
}
//______________________________________________________________________________
void NuVldMainFrame::PrintCurrentDBTable(void)
{
  NuVldUserData * user_data = NuVldUserData::Instance();

  GuiTablePrinter printer;
  printer.ScaleXSecWithEnergy( this->ScaleWithEnergy() );
  
  /* neutrino scattering data */

  if (fTabSql->GetCurrent() == 0) {

      if( ! user_data->CurrDBTableIsNull() )
                                  printer.PrintXSecTable( user_data->NuXSec() );
      else {
        fStatusBar -> SetText( "pointer to DBTable<T> is null", 1 );
        fLog       -> AddLine( "pointer to DBTable<T> is null"    );
      }

  } else

  /* electron scattering data */

  if (fTabSql->GetCurrent() == 1) {

      if( ! user_data->CurrDBTableIsNull() )
                              printer.PrintXSecTable( user_data->ElDiffXSec() );
      else {

        fStatusBar -> SetText( "pointer to DBTable<T> is null", 1 );
        fLog       -> AddLine( "pointer to DBTable<T> is null"    );
      }
   }//e
}
//______________________________________________________________________________
void NuVldMainFrame::RunFitter(void)
{
  LOG("NuVld", pDEBUG) << "Running fitter";

  // create a fit-kernel ans set GUI fit options

  fFitKernel->Reset();

  fFitKernel->SetFitParams(fNGFP);
  fFitKernel->SetFitRange(fXMinNmE->GetNumber(), fXMaxNmE->GetNumber());
  fFitKernel->SetScaleWithEnergy( this->ScaleWithEnergy() );

  fFitKernel->PrintConfig();

  NuVldUserData * user_data = NuVldUserData::Instance();

  // check whether a fitter was selected

  int entry_id = fFitterCBx->GetSelectedEntry()->EntryId();

  string fitter = kFitters[entry_id];

  if( fFitterCBx->GetSelectedEntry() && (strcmp(fitter.c_str(),"NONE") != 0) ) {

    // check the number of parameters to fit

    if( fNGFP->NFittedParams() > 0) {

      // if the fit is on stacked data & retrieve them and set as current

      if(fUseStackedChkB->GetState() == kButtonDown) {
         if(fTableStackCBx->GetSelectedEntry()) {

             int entry_id = fTableStackCBx->GetSelectedEntry()->EntryId();

             string xsec_table_name = fStackHandler->StackedDBTableName(entry_id);
             user_data->SetStackedDBTableAsCurr(xsec_table_name);

         } else new MsgBox(gClient->GetRoot(), fMain,  380, 250, kVerticalFrame,
                        " If you want to use a stacked data-set, then select one!");
      }

      // fit current the current data-set

      if( ! user_data->CurrDBTableIsNull() ) {

          LOG("NuVld", pDEBUG) << "Selected fitter: " << fitter;

          if      ( strcmp(fitter.c_str(),"SIMPLE") == 0 )     fFitKernel->DoSimpleFit();
          else if ( strcmp(fitter.c_str(),"NORM-FLOAT") == 0 ) fFitKernel->DoFloatingNormFit();
          else if ( strcmp(fitter.c_str(),"SIMPLE/UWF") == 0 )     {}
          else if ( strcmp(fitter.c_str(),"NORM-FLOAT/UWF") == 0 ) {}
          else {
            // should never be here
          }

          this->RunPostFitProcessor();

      } else new MsgBox(gClient->GetRoot(), fMain,
                 380, 250, kVerticalFrame, " You have selected no data to fit");

    } else new MsgBox(gClient->GetRoot(), fMain,
           380, 250, kVerticalFrame, " You have selected no parameters to fit");

  } else new MsgBox(gClient->GetRoot(), fMain,
                  380, 250, kVerticalFrame, " First, you must select a fitter");
}
//______________________________________________________________________________
void NuVldMainFrame::RunMcScanner(void)
{
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
    
  } else new MsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
             " You have not configured the multi-D MC param scanner properly "); 
}
//______________________________________________________________________________
void NuVldMainFrame::Run2dScanner(void)
{
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

  } else new MsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                           " You have not configured the 2d-scanner properly ");
}
//______________________________________________________________________________
void NuVldMainFrame::Run1dScanner(void)
{
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
/*
  double xmin = fFitKernel->chisq1d->GetX()[TMath::LocMin(fFitKernel->chisq1d->GetN(), fFitKernel->chisq1d->GetX())];
  double xmax = fFitKernel->chisq1d->GetX()[TMath::LocMax(fFitKernel->chisq1d->GetN(), fFitKernel->chisq1d->GetX())];
  double ymin = fFitKernel->chisq1d->GetY()[TMath::LocMin(fFitKernel->chisq1d->GetN(), fFitKernel->chisq1d->GetY())];
  double ymax = fFitKernel->chisq1d->GetY()[TMath::LocMax(fFitKernel->chisq1d->GetN(), fFitKernel->chisq1d->GetY())];

  TH1F * hframe = fFitTabFuncEmbCnv->GetCanvas()->DrawFrame(xmin,ymin,xmax,ymax);

  hframe->Draw();
*/
  fFitKernel->chisq1d->Draw("ALP");

  fLtxAuth->Draw();
  fLtxLink->Draw();
  
  fFitTabChisqEmbCnv->GetCanvas()->Update();
}
//______________________________________________________________________________
void NuVldMainFrame::PlotXSecBoundaries(TCanvas * c, bool clear)
{
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
}
//______________________________________________________________________________
void NuVldMainFrame::RunPostFitProcessor(void)
{
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
    renderer.SetErrorOption(this->ReadXSecSelectionListbox());
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
}
//______________________________________________________________________________
void NuVldMainFrame::PrintFitParameters(void)
{
  LOG("NuVld", pDEBUG) << "Getting fit func & and printing fit params";

  TF1 * func = fFitKernel -> FitFunction();

  if (func) {

    fFitTxtResults->AddLine("-----------------------------------------------");
  
    for(int i = 0; i < kNNGFitParams; i++) {
       if( fNGFP->IsFitted(i) ) {
         
          ostringstream fitparam;

          fitparam << fNGFP->ParamAsString(i) << " = "
                   << func->GetParameter(i) << " +/- " << func->GetParError(i);

          fFitTxtResults->AddLine(fitparam.str().c_str());
       }
    }

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
}
//______________________________________________________________________________
void NuVldMainFrame::DrawResiduals(void)
{
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
}
//______________________________________________________________________________
