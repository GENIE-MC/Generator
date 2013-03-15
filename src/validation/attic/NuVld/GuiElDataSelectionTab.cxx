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

#include <TGListBox.h>
#include <TGComboBox.h>
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
#include "ValidationTools/NuVld/GuiElDataSelectionTab.h"
#include "ValidationTools/NuVld/GuiElDataSelectionTabConstants.h"
#include "ValidationTools/NuVld/GuiSysLogSingleton.h"

using std::ostringstream;

using namespace genie;
using namespace genie::utils::str;
using namespace genie::nuvld;
using namespace genie::nuvld::constants;

ClassImp(GuiElDataSelectionTab)

//______________________________________________________________________________
GuiElDataSelectionTab::GuiElDataSelectionTab(DBConnection * db):
GuiDataSelectionDialog()
{
  fDBC = db;
}
//______________________________________________________________________________
GuiElDataSelectionTab::~GuiElDataSelectionTab()
{

}
//______________________________________________________________________________
TGCompositeFrame * GuiElDataSelectionTab::Create(
                                   TGCompositeFrame * tf, int width, int height)
{
  TGCompositeFrame * fTabElSql =
                        new TGCompositeFrame(tf, width, height, kVerticalFrame);

  fElExpGrpFrame = new TGGroupFrame(fTabElSql, "Experiment", kVerticalFrame);
  fElTgGrpFrm    = new TGGroupFrame(fTabElSql, "Target",     kVerticalFrame);

  fElExpLBx = new TGListBox(fElExpGrpFrame,  222);
  fElTgtLBx = new TGListBox(fElTgGrpFrm,     223);

  utils::gui::FillListBox( fElExpLBx,  kElExperiment );
  utils::gui::FillListBox( fElTgtLBx,  kElTarget     );

  fElExpLBx -> Resize (100,  60);
  fElTgtLBx -> Resize (100,  50);

  fElExpLBx -> SetMultipleSelections( true );
  fElTgtLBx -> SetMultipleSelections( true );

  fAllElExpChkB  = new TGCheckButton(fElExpGrpFrame, "Select all", 401);
  fAllElTgtChkB  = new TGCheckButton(fElTgGrpFrm,    "Select all", 402);

  fAllElExpChkB  -> Connect("Clicked()",
                   "genie::nuvld::GuiElDataSelectionTab", this,"SelectAllExp()");
  fAllElTgtChkB  -> Connect("Clicked()",
                   "genie::nuvld::GuiElDataSelectionTab", this,"SelectAllTargets()");

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

  utils::gui::FillComboBox( fElDrawXCBx, kElVarFrameName );

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
string GuiElDataSelectionTab::BundleSelectionsInString(void)
{
  ostringstream options;

  options << "KEY-LIST:" << this->BundleKeyListInString() << "$"
          << "CUTS:"     << this->BundleCutsInString()    << "$"
          << "DRAW_OPT:" << this->BundleDrawOptInString() << "$"
          << "DB-TYPE:eN-Diff-XSec";

  return options.str();
}
//______________________________________________________________________________
string GuiElDataSelectionTab::BundleKeyListInString(void)
{
  // Read experiment name selections
  string exprm = utils::gui::ListBoxSelectionAsString(fElExpLBx, kElExperiment);

  // Read target selections
  string targets = utils::gui::ListBoxSelectionAsString(fElTgtLBx, kElTarget);

  // Build key list
  string key_list = SqlUtils::build_e_key_list(fDBC->SqlServer(), exprm, targets);

  return key_list;
}
//______________________________________________________________________________
string GuiElDataSelectionTab::BundleCutsInString(void)
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
string GuiElDataSelectionTab::BundleDrawOptInString(void)
{
  int selected = fElDrawXCBx->GetSelected();

  const char * plot_var = kElVarMySQLName[selected];

  ostringstream draw_opt;

  draw_opt << "plot-var=" << plot_var;

  return draw_opt.str();
}
//______________________________________________________________________________
void GuiElDataSelectionTab::ResetSelections(void)
{
  utils::gui::ResetAllListBoxSelections( fElExpLBx );
  utils::gui::ResetAllListBoxSelections( fElTgtLBx );

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
void GuiElDataSelectionTab::SelectAllExp(void)
{
  if(fAllElExpChkB->GetState() == kButtonDown)
                                  utils::gui::SelectAllListBoxEntries(fElExpLBx);
  else utils::gui::ResetAllListBoxSelections(fElExpLBx);

  fElExpLBx->SelectionChanged();

  gClient->NeedRedraw(fElExpLBx->GetContainer());
}
//______________________________________________________________________________
void GuiElDataSelectionTab::SelectAllTargets(void)
{
  if(fAllElTgtChkB->GetState() == kButtonDown)
                                 utils::gui::SelectAllListBoxEntries(fElTgtLBx);
  else utils::gui::ResetAllListBoxSelections(fElTgtLBx);

  fElTgtLBx->SelectionChanged();

  gClient->NeedRedraw(fElTgtLBx->GetContainer());
}
//______________________________________________________________________________
