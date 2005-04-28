//_____________________________________________________________________________
/*!

\class    genie::nuvld::eDataSelectionTab

\brief    Electron Cross Section Data Selection Graphical Tab

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  November 01, 2004
*/
//_____________________________________________________________________________

#include <cassert>
#include <sstream>

#include <TGListBox.h>
#include <TGComboBox.h>
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
#include "NuVldGUI/eDataSelectionTab.h"
#include "NuVldGUI/eDataSelectionTabConstants.h"
#include "NuVldGUI/SysLogSingleton.h"
#include "Utils/GUIUtils.h"
#include "Utils/StringUtils.h"

using std::ostringstream;

using namespace genie;
using namespace genie::string_utils;
using namespace genie::nuvld;
using namespace genie::nuvld::constants;

ClassImp(eDataSelectionTab)

//______________________________________________________________________________
eDataSelectionTab::eDataSelectionTab(DBConnection * db):
DataSelectionDialog()
{  
  fDBC = db;
}
//______________________________________________________________________________
eDataSelectionTab::~eDataSelectionTab()
{

}
//______________________________________________________________________________
TGCompositeFrame * eDataSelectionTab::Create(
                                   TGCompositeFrame * tf, int width, int height)
{
  TGCompositeFrame * fTabElSql =
                        new TGCompositeFrame(tf, width, height, kVerticalFrame);

  fElExpGrpFrame = new TGGroupFrame(fTabElSql, "Experiment", kVerticalFrame);
  fElTgGrpFrm    = new TGGroupFrame(fTabElSql, "Target",     kVerticalFrame);

  fElExpLBx = new TGListBox(fElExpGrpFrame,  222);
  fElTgtLBx = new TGListBox(fElTgGrpFrm,     223);

  gui_utils::FillListBox( fElExpLBx,  kElExperiment );
  gui_utils::FillListBox( fElTgtLBx,  kElTarget     );

  fElExpLBx -> Resize (100,  60);
  fElTgtLBx -> Resize (100,  50);

  fElExpLBx -> SetMultipleSelections( true );
  fElTgtLBx -> SetMultipleSelections( true );

  fAllElExpChkB  = new TGCheckButton(fElExpGrpFrame, "Select all", 401);
  fAllElTgtChkB  = new TGCheckButton(fElTgGrpFrm,    "Select all", 402);

  fAllElExpChkB  -> Connect("Clicked()",
                   "genie::nuvld::eDataSelectionTab", this,"SelectAllExp()");
  fAllElTgtChkB  -> Connect("Clicked()",
                   "genie::nuvld::eDataSelectionTab", this,"SelectAllTargets()");

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

  this->SelectAllExp();
  this->SelectAllTargets();

  return fTabElSql;
}
//______________________________________________________________________________
string eDataSelectionTab::BundleSelectionsInString(void)
{
  ostringstream options;

  options << "KEY-LIST:" << this->BundleKeyListInString() << "$"
          << "CUTS:"     << this->BundleCutsInString()    << "$"
          << "DRAW_OPT:" << this->BundleDrawOptInString() << "$"
          << "DB-TYPE:eN-Diff-XSec";

  return options.str();
}
//______________________________________________________________________________
string eDataSelectionTab::BundleKeyListInString(void)
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
string eDataSelectionTab::BundleCutsInString(void)
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
string eDataSelectionTab::BundleDrawOptInString(void)
{
  int selected = fElDrawXCBx->GetSelected();

  const char * plot_var = kElVarMySQLName[selected];

  ostringstream draw_opt;

  draw_opt << "plot-var=" << plot_var;

  return draw_opt.str();
}
//______________________________________________________________________________
void eDataSelectionTab::ResetSelections(void)
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

  this->SelectAllExp();
  this->SelectAllTargets();
}
//______________________________________________________________________________
void eDataSelectionTab::SelectAllExp(void)
{
  if(fAllElExpChkB->GetState() == kButtonDown)
                                  gui_utils::SelectAllListBoxEntries(fElExpLBx);
  else gui_utils::ResetAllListBoxSelections(fElExpLBx);

  fElExpLBx->SelectionChanged();

  gClient->NeedRedraw(fElExpLBx->GetContainer());
}
//______________________________________________________________________________
void eDataSelectionTab::SelectAllTargets(void)
{
  if(fAllElTgtChkB->GetState() == kButtonDown)
                                 gui_utils::SelectAllListBoxEntries(fElTgtLBx);
  else gui_utils::ResetAllListBoxSelections(fElTgtLBx);

  fElTgtLBx->SelectionChanged();

  gClient->NeedRedraw(fElTgtLBx->GetContainer());
}
//______________________________________________________________________________
