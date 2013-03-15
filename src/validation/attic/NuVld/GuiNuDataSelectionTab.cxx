//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Nov 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

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

#include "Messenger/Messenger.h"
#include "Utils/GUIUtils.h"
#include "Utils/StringUtils.h"
#include "ValidationTools/NuVld/SqlUtils.hh"
#include "ValidationTools/NuVld/DBConnection.h"
#include "ValidationTools/NuVld/GuiNuDataSelectionTab.h"
#include "ValidationTools/NuVld/GuiNuDataSelectionTabConstants.h"
#include "ValidationTools/NuVld/GuiNuDataSelectionDialog.h"
#include "ValidationTools/NuVld/GuiNuMeasurementListDialog.h"
#include "ValidationTools/NuVld/GuiSysLogSingleton.h"
#include "ValidationTools/NuVld/GuiMsgBox.h"

using std::ostringstream;

using namespace genie;
using namespace genie::utils::str;
using namespace genie::nuvld;
using namespace genie::nuvld::constants;

ClassImp(GuiNuDataSelectionTab)

//______________________________________________________________________________
GuiNuDataSelectionTab::GuiNuDataSelectionTab(TGMainFrame * main, DBConnection * db):
GuiDataSelectionDialog()
{
  fMain = main;
  fDBC  = db;

  fPopupDialogLAM = false;
}
//______________________________________________________________________________
GuiNuDataSelectionTab::~GuiNuDataSelectionTab()
{

}
//______________________________________________________________________________
TGCompositeFrame * GuiNuDataSelectionTab::Create(
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

  utils::gui::FillListBox( fNuXSecErrLBx, kXSecErrType    );
  utils::gui::FillListBox( fNuExpLBx,     kExperimentName );
  utils::gui::FillListBox( fNuProcLBx,    kProcName       );
  utils::gui::FillListBox( fNuTypeLBx,    kNuType         );
  utils::gui::FillListBox( fNuTgtLBx,     kTarget         );

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
                      "genie::nuvld::GuiNuDataSelectionTab", this,"SelectAllExp()");
  fAllNuProcChkB -> Connect("Clicked()",
                     "genie::nuvld::GuiNuDataSelectionTab", this,"SelectAllXSec()");
  fAllNuTypesChkB-> Connect("Clicked()",
                   "genie::nuvld::GuiNuDataSelectionTab", this,"SelectAllProbes()");
  fAllNuTgtChkB  -> Connect("Clicked()",
                  "genie::nuvld::GuiNuDataSelectionTab", this,"SelectAllTargets()");

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
                   "genie::nuvld::GuiNuDataSelectionTab", this, "PopupNuGuiDataSelectionDialog()");
  fShowExpertNuDialogTBtn->Connect("Clicked()",
                 "genie::nuvld::GuiNuDataSelectionTab", this, "PopupNuXmlMeasurementListDialog()");

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
void GuiNuDataSelectionTab::SelectAllExp(void)
{
  if(fAllNuExpChkB->GetState() == kButtonDown)
                                  utils::gui::SelectAllListBoxEntries(fNuExpLBx);
  else utils::gui::ResetAllListBoxSelections(fNuExpLBx);

  fNuExpLBx->SelectionChanged();

  gClient->NeedRedraw(fNuExpLBx->GetContainer());
}
//______________________________________________________________________________
void GuiNuDataSelectionTab::SelectAllXSec(void)
{
  if(fAllNuProcChkB->GetState() == kButtonDown)
                                 utils::gui::SelectAllListBoxEntries(fNuProcLBx);
  else utils::gui::ResetAllListBoxSelections(fNuProcLBx);

  fNuProcLBx->SelectionChanged();

  gClient->NeedRedraw(fNuProcLBx->GetContainer());
}
//______________________________________________________________________________
void GuiNuDataSelectionTab::SelectAllProbes(void)
{
  if(fAllNuTypesChkB->GetState() == kButtonDown)
                                 utils::gui::SelectAllListBoxEntries(fNuTypeLBx);
  else utils::gui::ResetAllListBoxSelections(fNuTypeLBx);

  fNuTypeLBx->SelectionChanged();

  gClient->NeedRedraw(fNuTypeLBx->GetContainer());
}
//______________________________________________________________________________
void GuiNuDataSelectionTab::SelectAllTargets(void)
{
  if(fAllNuTgtChkB->GetState() == kButtonDown)
                                  utils::gui::SelectAllListBoxEntries(fNuTgtLBx);
  else utils::gui::ResetAllListBoxSelections(fNuTgtLBx);

  fNuTgtLBx->SelectionChanged();

  gClient->NeedRedraw(fNuTgtLBx->GetContainer());
}
//______________________________________________________________________________
string GuiNuDataSelectionTab::BundleSelectionsInString(void)
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
string GuiNuDataSelectionTab::BundleKeyListInString(void)
{
  if(fPopupDialogLAM) return fPopupDialog->BundleKeyListInString();

  bool is_connected;
  if( !fDBC->SqlServer() ) is_connected = false;
  else {
    is_connected = fDBC->SqlServer()->IsConnected();
  }
  if(!is_connected) return "";

  // Read experiment name selections
  string experiments = utils::gui::ListBoxSelectionAsString(
                                                fNuExpLBx, kExperimentMySQLName);
  // Read xsec selections
  string xsecs = utils::gui::ListBoxSelectionAsString(fNuProcLBx, kProcMySQLName);
  // Read neutrino selections
  string nus = utils::gui::ListBoxSelectionAsString(fNuTypeLBx, kNuTypeMySQLName);
  // Read target selections
  string targets = utils::gui::ListBoxSelectionAsString(
                                                    fNuTgtLBx, kTargetMySQLName);

  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

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
string GuiNuDataSelectionTab::BundleCutsInString(void)
{
  if(fPopupDialogLAM) return fPopupDialog->BundleCutsInString();

  float Emin = fEMinNmE->GetNumber();
  float Emax = fEMaxNmE->GetNumber();

  ostringstream cuts;

  cuts << "Emin=" << Emin << ";" << "Emax=" << Emax;

  return cuts.str();
}
//______________________________________________________________________________
string GuiNuDataSelectionTab::BundleDrawOptInString(void)
{
  if(fPopupDialogLAM) return fPopupDialog->BundleDrawOptInString();

  bool   scale_e = (fScaleWithEvChkB->GetState() == kButtonDown);
  string err_opt = this->ReadXSecErrorListbox();

  ostringstream dopt;

  dopt << "scale-with-energy=";
  if(scale_e) dopt << "yes;";
  else        dopt << "no;";

  dopt << "err-opt=" << err_opt << ";";

  return dopt.str();
}
//______________________________________________________________________________
void GuiNuDataSelectionTab::ResetSelections(void)
{
  utils::gui::ResetAllListBoxSelections( fNuExpLBx  );
  utils::gui::ResetAllListBoxSelections( fNuProcLBx );
  utils::gui::ResetAllListBoxSelections( fNuTypeLBx );
  utils::gui::ResetAllListBoxSelections( fNuTgtLBx  );

  fEMinNmE->SetNumber(kEmin);
  fEMaxNmE->SetNumber(kEmax);

  fNuXSecErrLBx  -> Select (1);

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
string GuiNuDataSelectionTab::ReadXSecErrorListbox(void)
{
  ostringstream err;

  TGLBEntry * selected_entry = fNuXSecErrLBx->GetSelectedEntry();

  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

  if(selected_entry) {
    int sid = selected_entry->EntryId();
    err << kXSecErrType[sid] << "-noE";
    syslog->Log()->AddLine( Concat("XSec Errors Selection: ", sid) );
  } else {
    err << "allXsec-noE";
    syslog->Log()->AddLine("No Cross Section Error Selection - setting default");
  }
  LOG("NuVld", pDEBUG) << "error selection = " << err.str().c_str();

  return err.str();
}
//______________________________________________________________________________
void GuiNuDataSelectionTab::PopupNuGuiDataSelectionDialog(void)
{
  bool IsConnected;

  if( !fDBC->SqlServer() ) IsConnected = false;
  else IsConnected = fDBC->SqlServer()->IsConnected();

  if(IsConnected) {

     if(!fPopupDialogLAM) {

      fPopupDialog = new GuiNuDataSelectionDialog(
                  gClient->GetRoot(), fMain, &fPopupDialogLAM,
                                            750, 500, kHorizontalFrame, fDBC );
     } else {

       new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
           "Another selection dialog has locked my attention. Close it first.");
     }

  } else {
      new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                   "You must be connected to the data-base to use this option");
  }
}
//______________________________________________________________________________
void GuiNuDataSelectionTab::PopupNuXmlMeasurementListDialog(void)
{
  bool IsConnected;

  if( !fDBC->SqlServer() ) IsConnected = false;
  else IsConnected = fDBC->SqlServer()->IsConnected();

  if(IsConnected) {

     if(!fPopupDialogLAM) {

      fPopupDialog = new GuiNuMeasurementListDialog(
                  gClient->GetRoot(), fMain, &fPopupDialogLAM,
                                               650, 400, kVerticalFrame, fDBC );
     } else {

       new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
           "Another selection dialog has locked my attention. Close it first.");
     }

  } else {
      new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                   "You must be connected to the data-base to use this option");
  }
}
//______________________________________________________________________________
