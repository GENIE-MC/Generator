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

#include <TGFrame.h>
#include <TGListBox.h>
#include <TGComboBox.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGText.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSystem.h>

#include "Messenger/Messenger.h"
#include "Utils/GUIUtils.h"
#include "Utils/StringUtils.h"
#include "ValidationTools/NuVld/SqlUtils.hh"
#include "ValidationTools/NuVld/DBConnection.h"
#include "ValidationTools/NuVld/GuiSFDataSelectionTab.h"
#include "ValidationTools/NuVld/GuiSFDataSelectionTabConstants.h"
#include "ValidationTools/NuVld/GuiSysLogSingleton.h"

using std::ostringstream;

using namespace genie;
using namespace genie::utils::str;
using namespace genie::nuvld;
using namespace genie::nuvld::constants;

ClassImp(GuiSFDataSelectionTab)

//______________________________________________________________________________
GuiSFDataSelectionTab::GuiSFDataSelectionTab(DBConnection * db):
GuiDataSelectionDialog()
{
  fDBC = db;
}
//______________________________________________________________________________
GuiSFDataSelectionTab::~GuiSFDataSelectionTab()
{

}
//______________________________________________________________________________
TGCompositeFrame * GuiSFDataSelectionTab::Create(
                                   TGCompositeFrame * tf, int width, int height)
{
  UInt_t kv = kVerticalFrame;

  fTabSFSql = new TGCompositeFrame(tf, width, height, kv);

  fErrGrpFrm       = new TGGroupFrame(fTabSFSql, "Err Type",           kv);
  fExpGrpFrm       = new TGGroupFrame(fTabSFSql, "Experiment",         kv);
  fSFGrpFrm        = new TGGroupFrame(fTabSFSql, "SF, R=sT/sL",        kv);
  fKineGrpFrm      = new TGGroupFrame(fTabSFSql, "Kine: Q2,x min/max", kv);
  fInitStateGrpFrm = new TGGroupFrame(fTabSFSql, "Initial State",      kv);
  fPlotVarGrpFrm   = new TGGroupFrame(fTabSFSql, "Plot Variable:",     kv);

  fErrLBx     = new TGListBox  (fErrGrpFrm,         2);
  fExpLBx     = new TGListBox  (fExpGrpFrm,         2);
  fSFLBx      = new TGListBox  (fSFGrpFrm,          2);
  fRLBx       = new TGListBox  (fSFGrpFrm,          2);
  fProbeLBx   = new TGListBox  (fInitStateGrpFrm,   2);
  fTgtLBx     = new TGListBox  (fInitStateGrpFrm,   2);
  fPlotVarCBx = new TGComboBox (fPlotVarGrpFrm,     2);
  //fSFxLBx   = new TGListBox(fKineGrpFrm, 2);

  utils::gui::FillListBox  ( fErrLBx,     kSFErrType        );
  utils::gui::FillListBox  ( fExpLBx,     kSFExperimentName );
  utils::gui::FillListBox  ( fSFLBx,      kSFName           );
  utils::gui::FillListBox  ( fRLBx,       kSFR              );
  utils::gui::FillListBox  ( fProbeLBx,   kSFProbe          );
  utils::gui::FillListBox  ( fTgtLBx,     kSFTarget         );
  utils::gui::FillComboBox ( fPlotVarCBx, kSFPlotVar        );

  fErrLBx     -> Resize (100,  60);
  fExpLBx     -> Resize (100,  60);
  fSFLBx      -> Resize (100,  40);
  fProbeLBx   -> Resize (100,  60);
  fTgtLBx     -> Resize (100,  60);
  fRLBx       -> Resize (100,  60);
  fPlotVarCBx -> Resize (100,  20);
  //fSFxLBx   -> Resize (100,  60);

  fErrLBx   -> SetMultipleSelections( false );
  fExpLBx   -> SetMultipleSelections( true  );
  fSFLBx    -> SetMultipleSelections( false );
  fProbeLBx -> SetMultipleSelections( true  );
  fTgtLBx   -> SetMultipleSelections( true  );
  fRLBx     -> SetMultipleSelections( true  );
  //fSFxLBx   -> SetMultipleSelections( true  );

  fAllExpChkB    = new TGCheckButton(fExpGrpFrm,       "Select all", 371);
  fAllProbesChkB = new TGCheckButton(fInitStateGrpFrm, "Select all", 372);
  fAllTgtChkB    = new TGCheckButton(fInitStateGrpFrm, "Select all", 373);

  fAllExpChkB    -> Connect("Clicked()",
                     "genie::nuvld::GuiSFDataSelectionTab", this,"SelectAllExp()");
  fAllProbesChkB -> Connect("Clicked()",
                     "genie::nuvld::GuiSFDataSelectionTab",this,"SelectAllProbes()");
  fAllTgtChkB    -> Connect("Clicked()",
                     "genie::nuvld::GuiSFDataSelectionTab",this,"SelectAllTargets()");

  // kinematical variables: x,Q2 range

  TGNumberFormat::EStyle style = TGNumberFormat::kNESReal;

  fMinQ2NmE  = new TGNumberEntry(fKineGrpFrm,   0,  12, 3, style);
  fMaxQ2NmE  = new TGNumberEntry(fKineGrpFrm, 100., 12, 3, style);
  fMinXNmE   = new TGNumberEntry(fKineGrpFrm,   0., 12, 5, style);
  fMaxXNmE   = new TGNumberEntry(fKineGrpFrm,   1., 12, 5, style);

  //fSFLoadxTBtn = new TGTextButton (fKineGrpFrm, "Load x... ", 374);
  //fSFLoadxTBtn->Connect("Clicked()",
  //                    "genie::nuvld::GuiSFDataSelectionTab", this, "SFLoadx()");

  fErrGrpFrm       -> AddFrame ( fErrLBx        );
  fExpGrpFrm       -> AddFrame ( fExpLBx        );
  fExpGrpFrm       -> AddFrame ( fAllExpChkB    );
  fSFGrpFrm        -> AddFrame ( fSFLBx         );
  fSFGrpFrm        -> AddFrame ( fRLBx          );
  fInitStateGrpFrm -> AddFrame ( fProbeLBx      );
  fInitStateGrpFrm -> AddFrame ( fAllProbesChkB );
  fInitStateGrpFrm -> AddFrame ( fTgtLBx        );
  fInitStateGrpFrm -> AddFrame ( fAllTgtChkB    );
  fKineGrpFrm      -> AddFrame ( fMinQ2NmE      );
  fKineGrpFrm      -> AddFrame ( fMaxQ2NmE      );
  fKineGrpFrm      -> AddFrame ( fMinXNmE       );
  fKineGrpFrm      -> AddFrame ( fMaxXNmE       );
  fPlotVarGrpFrm   -> AddFrame ( fPlotVarCBx    );
  //fKineGrpFrm      -> AddFrame ( fSFxLBx          );
  //fKineGrpFrm      -> AddFrame ( fSFLoadxTBtn     );

  fTabSFSql -> AddFrame( fErrGrpFrm        );
  fTabSFSql -> AddFrame( fExpGrpFrm        );
  fTabSFSql -> AddFrame( fSFGrpFrm         );
  fTabSFSql -> AddFrame( fInitStateGrpFrm  );
  fTabSFSql -> AddFrame( fKineGrpFrm       );
  fTabSFSql -> AddFrame( fPlotVarGrpFrm    );

  this->ResetSelections();

  return fTabSFSql;
}
//______________________________________________________________________________
void GuiSFDataSelectionTab::SelectAllExp(void)
{
  if(fAllExpChkB->GetState() == kButtonDown)
                                  utils::gui::SelectAllListBoxEntries(fExpLBx);
  else utils::gui::ResetAllListBoxSelections(fExpLBx);

  fExpLBx->SelectionChanged();

  gClient->NeedRedraw(fExpLBx->GetContainer());
}
//______________________________________________________________________________
void GuiSFDataSelectionTab::SelectAllProbes(void)
{
  if(fAllProbesChkB->GetState() == kButtonDown)
                                 utils::gui::SelectAllListBoxEntries(fProbeLBx);
  else utils::gui::ResetAllListBoxSelections(fProbeLBx);

  fProbeLBx->SelectionChanged();

  gClient->NeedRedraw(fProbeLBx->GetContainer());
}
//______________________________________________________________________________
void GuiSFDataSelectionTab::SelectAllTargets(void)
{
  if(fAllTgtChkB->GetState() == kButtonDown)
                                  utils::gui::SelectAllListBoxEntries(fTgtLBx);
  else utils::gui::ResetAllListBoxSelections(fTgtLBx);

  fTgtLBx->SelectionChanged();

  gClient->NeedRedraw(fTgtLBx->GetContainer());
}
//______________________________________________________________________________
void GuiSFDataSelectionTab::ResetSelections(void)
{
  utils::gui::ResetAllListBoxSelections( fExpLBx   );
  utils::gui::ResetAllListBoxSelections( fSFLBx      );
  utils::gui::ResetAllListBoxSelections( fRLBx     );
  utils::gui::ResetAllListBoxSelections( fProbeLBx );
  utils::gui::ResetAllListBoxSelections( fTgtLBx   );

  fMinQ2NmE      -> SetNumber (0);
  fMaxQ2NmE      -> SetNumber (100);
  fMinXNmE       -> SetNumber (0);
  fMaxXNmE       -> SetNumber (1);
  fErrLBx        -> Select    (2);
  fSFLBx         -> Select    (0);
  fRLBx          -> Select    (0);
  fRLBx          -> Select    (2);
  fPlotVarCBx    -> Select    (0);
  fAllExpChkB    -> SetOn     (kTRUE);
  fAllProbesChkB -> SetOn     (kTRUE);
  fAllTgtChkB    -> SetOn     (kTRUE);

  this -> SelectAllExp();
  this -> SelectAllProbes();
  this -> SelectAllTargets();
}
//______________________________________________________________________________
string GuiSFDataSelectionTab::BundleSelectionsInString(void)
{
  ostringstream options;

  options << "KEY-LIST:" << this->BundleKeyListInString()  << "$"
          << "CUTS:"     << this->BundleCutsInString()     << "$"
          << "DRAW_OPT:" << this->BundleDrawOptInString()  << "$"
          << "DB-TYPE:SF";

  return options.str();
}
//______________________________________________________________________________
string GuiSFDataSelectionTab::BundleKeyListInString(void)
{
  // Read experiment name selections
  string experiments = utils::gui::ListBoxSelectionAsString(
                                                  fExpLBx, kSFExperimentName);
  // Read SF selection
  string sf = utils::gui::ListBoxSelectionAsString(fSFLBx, kSFName);
  // Read probe selections
  string probes = utils::gui::ListBoxSelectionAsString(fProbeLBx, kSFProbe);
  // Read target selections
  string targets = utils::gui::ListBoxSelectionAsString(fTgtLBx, kSFTarget);

  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

  syslog->Log()->AddLine( "Requested key : ");
  syslog->Log()->AddLine( Concat("- Experiments... : ", experiments.c_str()) );
  syslog->Log()->AddLine( Concat("- SF............ : ", sf.c_str())          );
  syslog->Log()->AddLine( Concat("- Probes........ : ", probes.c_str())      );
  syslog->Log()->AddLine( Concat("- Targets....... : ", targets.c_str())     );

  // Build key list
  string key_list = SqlUtils::build_sf_key_list(
                          fDBC->SqlServer(), experiments, sf, probes, targets);

  return key_list;
}
//______________________________________________________________________________
string GuiSFDataSelectionTab::BundleCutsInString(void)
{
  float Q2min = fMinQ2NmE -> GetNumber();
  float Q2max = fMaxQ2NmE -> GetNumber();
  float xmin  = fMinXNmE  -> GetNumber();
  float xmax  = fMaxXNmE  -> GetNumber();

  // Read R selection
  string R = utils::gui::ListBoxSelectionAsString(fRLBx, kSFR);

  // Read x selections
  //string x = utils::gui::ListBoxSelectionAsString(fSFxLBx, kSFR);

  ostringstream cuts;

  cuts << "Q2min=" << Q2min << ";" << "Q2max=" << Q2max << ";";
  cuts << "xmin="  << xmin  << ";" << "xmax="  << xmax  << ";";
  cuts << "R=" << R;

  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

  syslog->Log()->AddLine( Concat("Requested cuts : ", cuts.str().c_str()) );

  return cuts.str();
}
//______________________________________________________________________________
string GuiSFDataSelectionTab::BundleDrawOptInString(void)
{
  int selected = fPlotVarCBx->GetSelected();

  const char * plot_var = kSFPlotVar[selected];

  ostringstream draw_opt;

  draw_opt << "plot-var=" << plot_var;

  return draw_opt.str();
}
//______________________________________________________________________________
/*
void GuiSFDataSelectionTab::SFLoadx(void)
{
  bool IsConnected;

  if( !fDBC->SqlServer() ) IsConnected = false;
  else IsConnected = fDBC->SqlServer()->IsConnected();

  if(IsConnected) {

     string query = "SELECT DISTINCT x from STRUCTURE_FUNCTION";
     TSQLResult * res = fDBC->SqlServer()->Query(query.c_str());

     const int nrows = res->GetRowCount();
     vector<string> x(nrows);

     for (int i = 0; i < nrows; i++) {

       TSQLRow * row = res->Next();
       x[i] = row->GetField(0);

       LOG("NuVld", pINFO)
              << "Adding x in SF kinematics: " << "x[" << i << "] = " << x[i];
       delete row;
     }
     delete res;

     utils::gui::FillListBox(fSFxLBx,&x);

     fSFxLBx->MapSubwindows();
     fSFxLBx->Layout();

     fSFxLBx->SelectionChanged();
     gClient->NeedRedraw(fSFxLBx->GetContainer());

     gSystem->ProcessEvents();
  }
}
*/
//______________________________________________________________________________
