//_____________________________________________________________________________
/*!

\class    genie::nuvld::vDataSelectionTab

\brief    Neutrino Cross Section Data Selection Graphical Tab

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  November 01, 2004
*/
//_____________________________________________________________________________

#include <cassert>
#include <sstream>

#include <TGClient.h>
#include <TGFrame.h>
#include <TGListBox.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TGText.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include "DBUtils/SqlUtils.hh"
#include "Messenger/Messenger.h"
#include "NuVldGUI/DBConnection.h"
#include "NuVldGUI/vDataSelectionTab.h"
#include "NuVldGUI/vDataSelectionTabConstants.h"
#include "NuVldGUI/vDataSelectionDialog.h"
#include "NuVldGUI/vMeasurementListDialog.h"
#include "NuVldGUI/SysLogSingleton.h"
#include "NuVldGUI/MsgBox.h"
#include "Utils/GUIUtils.h"
#include "Utils/StringUtils.h"

using std::ostringstream;

using namespace genie;
using namespace genie::string_utils;
using namespace genie::nuvld;
using namespace genie::nuvld::constants;

ClassImp(vDataSelectionTab)

//______________________________________________________________________________
vDataSelectionTab::vDataSelectionTab(TGMainFrame * main, DBConnection * db):
DataSelectionDialog()
{
  fMain = main;
  fDBC  = db;
  
  fPopupDialogLAM = false;
}
//______________________________________________________________________________
vDataSelectionTab::~vDataSelectionTab()
{

}
//______________________________________________________________________________
TGCompositeFrame * vDataSelectionTab::Create(
                                   TGCompositeFrame * tf, int width, int height)
{
  UInt_t kv = kVerticalFrame;
  
  TGCompositeFrame * fTabNuSql = new TGCompositeFrame(tf, width, height, kv);

  fNuXSecErrGrpFrm   = new TGGroupFrame(fTabNuSql, "Cross Section Err",  kv);
  fNuExpGrpFrm       = new TGGroupFrame(fTabNuSql, "Experiment",         kv);
  fNuXSecGrpFrm      = new TGGroupFrame(fTabNuSql, "Cross Section",      kv);
  fEnergyGrpFrm      = new TGGroupFrame(fTabNuSql, "Energy Range (GeV)", kv);
  fNuInitStateGrpFrm = new TGGroupFrame(fTabNuSql, "Initial State",      kv);

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

  fAllNuExpChkB  -> Connect("Clicked()",
                      "genie::nuvld::vDataSelectionTab", this,"SelectAllExp()");
  fAllNuProcChkB -> Connect("Clicked()",
                     "genie::nuvld::vDataSelectionTab", this,"SelectAllXSec()");
  fAllNuTypesChkB-> Connect("Clicked()",
                   "genie::nuvld::vDataSelectionTab", this,"SelectAllProbes()");
  fAllNuTgtChkB  -> Connect("Clicked()",
                  "genie::nuvld::vDataSelectionTab", this,"SelectAllTargets()");

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

  TGNumberFormat::EStyle rstyle = TGNumberFormat::kNESReal;
  
  fEMinNmE = new TGNumberEntry(fEnergyGrpFrm, kEmin, 6, 1, rstyle);
  fEMaxNmE = new TGNumberEntry(fEnergyGrpFrm, kEmax, 6, 1, rstyle);

  fMinELb = new TGLabel(fEnergyGrpFrm, new TGString( "min:"));
  fMaxELb = new TGLabel(fEnergyGrpFrm, new TGString( "max:"));

  fEnergyGrpFrm -> AddFrame ( fMinELb  );
  fEnergyGrpFrm -> AddFrame ( fEMinNmE );
  fEnergyGrpFrm -> AddFrame ( fMaxELb  );
  fEnergyGrpFrm -> AddFrame ( fEMaxNmE );

  fNuTabBtnSpacerLb = new TGLabel(fTabNuSql, new TGString(" "));

  fShowFullNuDialogTBtn   = new TGTextButton (fTabNuSql, "More data selections... ", 76);
  fShowExpertNuDialogTBtn = new TGTextButton (fTabNuSql, "Expert mode...          ", 77);

  fShowFullNuDialogTBtn->Connect("Clicked()", 
                   "genie::nuvld::vDataSelectionTab", this, "PopupNuDataSelectionDialog()");
  fShowExpertNuDialogTBtn->Connect("Clicked()", 
                 "genie::nuvld::vDataSelectionTab", this, "PopupNuMeasurementListDialog()");

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

  return fTabNuSql;
}
//______________________________________________________________________________
void vDataSelectionTab::SelectAllExp(void)
{
  if(fAllNuExpChkB->GetState() == kButtonDown)
                                  gui_utils::SelectAllListBoxEntries(fNuExpLBx);
  else gui_utils::ResetAllListBoxSelections(fNuExpLBx);

  fNuExpLBx->SelectionChanged();

  gClient->NeedRedraw(fNuExpLBx->GetContainer());
}
//______________________________________________________________________________
void vDataSelectionTab::SelectAllXSec(void)
{
  if(fAllNuProcChkB->GetState() == kButtonDown)
                                 gui_utils::SelectAllListBoxEntries(fNuProcLBx);
  else gui_utils::ResetAllListBoxSelections(fNuProcLBx);

  fNuProcLBx->SelectionChanged();

  gClient->NeedRedraw(fNuProcLBx->GetContainer());
}
//______________________________________________________________________________
void vDataSelectionTab::SelectAllProbes(void)
{
  if(fAllNuTypesChkB->GetState() == kButtonDown)
                                 gui_utils::SelectAllListBoxEntries(fNuTypeLBx);
  else gui_utils::ResetAllListBoxSelections(fNuTypeLBx);

  fNuTypeLBx->SelectionChanged();

  gClient->NeedRedraw(fNuTypeLBx->GetContainer());
}
//______________________________________________________________________________
void vDataSelectionTab::SelectAllTargets(void)
{
  if(fAllNuTgtChkB->GetState() == kButtonDown)
                                  gui_utils::SelectAllListBoxEntries(fNuTgtLBx);
  else gui_utils::ResetAllListBoxSelections(fNuTgtLBx);

  fNuTgtLBx->SelectionChanged();

  gClient->NeedRedraw(fNuTgtLBx->GetContainer());
}
//______________________________________________________________________________
string vDataSelectionTab::BundleSelectionsInString(void)
{
  if(fPopupDialogLAM) {
    LOG("NuVld", pDEBUG) << "Found LAM flag from a popup dialog";

    return fPopupDialog->BundleSelectionsInString();
  }
  
  ostringstream options;

  options << "KEY-LIST:" << this->BundleKeyListInString()  << "$"
          << "CUTS:"     << this->BundleCutsInString()     << "$"
          << "DRAW_OPT:" << this->BundleDrawOptInString()  << "$"
          << "DB-TYPE:vN-XSec";

  return options.str();
}
//______________________________________________________________________________
string vDataSelectionTab::BundleKeyListInString(void)
{
  if(fPopupDialogLAM) return fPopupDialog->BundleKeyListInString();

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

  SysLogSingleton * syslog = SysLogSingleton::Instance();

  syslog->Log()->AddLine( Concat("requested experiments : ", experiments.c_str()) );
  syslog->Log()->AddLine( Concat("requested measurements : ", xsecs.c_str()) );
  syslog->Log()->AddLine( Concat("requested neutrino beams : ", nus.c_str()) );
  syslog->Log()->AddLine( Concat("requested targets : ", targets.c_str()) );

  // Build key list
  string key_list = SqlUtils::build_v_key_list(
                fDBC->SqlServer(), experiments, xsecs, nus, targets);

  return key_list;
}
//______________________________________________________________________________
string vDataSelectionTab::BundleCutsInString(void)
{
  if(fPopupDialogLAM) return fPopupDialog->BundleCutsInString();

  float Emin = fEMinNmE->GetNumber();
  float Emax = fEMaxNmE->GetNumber();

  ostringstream cuts;

  cuts << "Emin=" << Emin << ";" << "Emax=" << Emax;

  return cuts.str();
}
//______________________________________________________________________________
string vDataSelectionTab::BundleDrawOptInString(void)
{
  if(fPopupDialogLAM) return fPopupDialog->BundleDrawOptInString();

  if(fScaleWithEvChkB->GetState() == kButtonDown) return "scale-with-energy";
  else return "";
}
//______________________________________________________________________________
void vDataSelectionTab::ResetSelections(void)
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

  this->SelectAllExp();
  this->SelectAllXSec();
  this->SelectAllProbes();
  this->SelectAllTargets();
}
//______________________________________________________________________________
string vDataSelectionTab::ReadXSecSelectionListbox(void)
{
  ostringstream err;

  TGLBEntry * selected_entry = fNuXSecErrLBx->GetSelectedEntry();

  SysLogSingleton * syslog = SysLogSingleton::Instance();

  if(selected_entry) {

    err << kXSecErrDrawOpt[ selected_entry->EntryId() ] << "-noE";

      syslog->Log()->AddLine( Concat("XSec Errors Selection: ",
                                                  selected_entry->EntryId()) );
  } else {
    err << "allXsec-noE";
      syslog->Log()->AddLine(
                         "No Cross Section Error Selection - setting default" );
  }
  LOG("NuVld", pDEBUG) << "error selection = " << err.str().c_str();

  return err.str();
}
//______________________________________________________________________________
void vDataSelectionTab::PopupNuDataSelectionDialog(void)
{
  bool IsConnected;

  if( !fDBC->SqlServer() ) IsConnected = false;
  else IsConnected = fDBC->SqlServer()->IsConnected();

  if(IsConnected) {

     if(!fPopupDialogLAM) {

      fPopupDialog = new vDataSelectionDialog(
                  gClient->GetRoot(), fMain, &fPopupDialogLAM,
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
void vDataSelectionTab::PopupNuMeasurementListDialog(void)
{
  bool IsConnected;

  if( !fDBC->SqlServer() ) IsConnected = false;
  else IsConnected = fDBC->SqlServer()->IsConnected();

  if(IsConnected) {

     if(!fPopupDialogLAM) {

      fPopupDialog = new vMeasurementListDialog(
                  gClient->GetRoot(), fMain, &fPopupDialogLAM,
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
