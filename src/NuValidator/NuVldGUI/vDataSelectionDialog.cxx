//_____________________________________________________________________________
/*!

\class    genie::nuvld::vDataSelectionDialog

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  November 01, 2004
*/
//_____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>

#include <TGText.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include "NuVldGUI/vDataSelectionDialog.h"
#include "XmlParser/ParserUtils.h"
#include "Utils/GUIUtils.h"

using std::cout;
using std::endl;
using std::map;
using std::ostringstream;
using std::setw;
using std::setfill;

using namespace genie;
using namespace genie::nuvld;

ClassImp(vDataSelectionDialog)

const char * k_xsec_err_types[] = {"no error", "stat. only", "stat.+syst.", 0};

//______________________________________________________________________________
vDataSelectionDialog::vDataSelectionDialog(
             const TGWindow *p, const TGWindow *main, bool & attn,
                        UInt_t w, UInt_t h, UInt_t options, DBConnection * db):
DataSelectionDialog(attn)
{
  _db = db;
  
  _main = new TGTransientFrame(p, main, w, h, options);
  _main->Connect("CloseWindow()",
                  "genie::nuvld::vDataSelectionDialog", this, "CloseWindow()");

  _main_left_frame   = new TGCompositeFrame(_main, 3, 3, kVerticalFrame);
  _main_right_frame  = new TGCompositeFrame(_main, 3, 3, kVerticalFrame);
  
  ULong_t hintMLeftFrameLayout    = kLHintsCenterY;
  ULong_t hintMRightFrameLayout   = kLHintsTop | kLHintsExpandX | kLHintsExpandY;
  
  _mleft_frame_layout  = new TGLayoutHints(hintMLeftFrameLayout,    1, 1,  1, 1);
  _mright_frame_layout = new TGLayoutHints(hintMRightFrameLayout,   1, 1,  1, 1);

  this->BuildLeftFrameWidgets();
  this->BuildRightFrameWidgets();
        
  _main  -> AddFrame ( _main_left_frame,    _mleft_frame_layout  );
  _main  -> AddFrame ( _main_right_frame,   _mright_frame_layout );
  
  _main->MapSubwindows();
  _main->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  _main->SetWindowName("Neutrino Data Selection Dialog");

  _main->MapWindow();

  // init selection
  SelectAllExp();

  //gClient->WaitFor(_main);
}
//______________________________________________________________________________
vDataSelectionDialog::~vDataSelectionDialog()
{
  delete _mleft_frame_layout;
  delete _mright_frame_layout;
  delete _listbox_layout;
  delete _button_layout;
  delete _close_button;
  delete _W_cut;
  delete _E_min;
  delete _E_max;
  delete _E_minLabel;
  delete _E_maxLabel;
  delete _select_all_exp;
  delete _select_all_xsec;
  delete _select_all_nu;
  delete _select_all_target;
  delete _scale_with_energy;
  delete _measurements_listbox;
  delete _xsec_err_listbox;
  delete _exp_listbox;
  delete _xsec_listbox;
  delete _nu_listbox;
  delete _target_listbox;
  delete _xsec_err_group_frame; 
  delete _exp_group_frame;
  delete _obs_group_frame;
  //delete _energy_matrix_layout;
  delete _energy_group_frame;
  delete _wcut_group_frame;
  delete _cuts_group_frame;  
  delete _init_state_group_frame;
  delete _target_group_frame;
  delete _reaction_group_frame;  
  delete _main_left_frame;
  delete _main_right_frame;
  delete _main;
}
//______________________________________________________________________________
void vDataSelectionDialog::BuildLeftFrameWidgets(void)
{
  _xsec_err_group_frame   = new TGGroupFrame(_main_left_frame,
                                          "Cross Section Err",  kVerticalFrame);
  _exp_group_frame        = new TGGroupFrame(_main_left_frame,
                                          "Experiment",         kVerticalFrame);
  _cuts_group_frame       = new TGGroupFrame(_main_left_frame,
                                                        "Cuts", kVerticalFrame);
  _energy_group_frame     = new TGGroupFrame(_cuts_group_frame,
                                          "Energy Range (GeV)", kVerticalFrame);
  _wcut_group_frame       = new TGGroupFrame(_cuts_group_frame,
                                          "Wmax (GeV) for SPP", kVerticalFrame);

  _energy_matrix_layout = new TGMatrixLayout(_energy_group_frame, 0, 2, 2);
  _energy_group_frame->SetLayoutManager( _energy_matrix_layout );

  _xsec_err_listbox = new TGListBox(_xsec_err_group_frame, 2);
  _exp_listbox      = new TGListBox(_exp_group_frame,      2);

  this->LoadExperimentsFromDB();

  gui_utils::FillListBox( _xsec_err_listbox, k_xsec_err_types );
  
  _xsec_err_listbox -> Resize (120,  60);
  _exp_listbox      -> Resize (120,  90);

  _xsec_err_listbox -> SetMultipleSelections( false );
  _exp_listbox      -> SetMultipleSelections( true  );

  _select_all_exp    = new TGCheckButton(_exp_group_frame,     "Select all", 71);

  //--- FORCE CONSISTENCY BETWEEN SELECTED VALUES
  
  // Update widgets / force consistency every time the expts listbox is clicked
  
  _exp_listbox ->Connect("SelectionChanged()","genie::nuvld::vDataSelectionDialog",
                                       this,"MakeConsistentWithExpListbox()");
    
  //--- "Select all" action

  _select_all_exp    -> Connect("Clicked()","genie::nuvld::vDataSelectionDialog",
                                                       this,"SelectAllExp()");


  _E_min = new TGNumberEntry(
                   _energy_group_frame,   0.1, 6, 1, TGNumberFormat::kNESReal);
  _E_max = new TGNumberEntry(
                   _energy_group_frame, 200.0, 6, 1, TGNumberFormat::kNESReal);

  _E_minLabel = new TGLabel(_energy_group_frame, new TGString( "min:"));
  _E_maxLabel = new TGLabel(_energy_group_frame, new TGString( "max:"));


  _W_cut = new TGNumberEntry(
                     _wcut_group_frame,   1.4, 6, 1, TGNumberFormat::kNESReal);

                   
  _scale_with_energy  = new TGCheckButton(
                                      _main_left_frame, "Scale With Energy", 75);
    
  _xsec_err_group_frame  -> AddFrame( _xsec_err_listbox  );
  _exp_group_frame       -> AddFrame( _exp_listbox       );
  _exp_group_frame       -> AddFrame( _select_all_exp    );

  _energy_group_frame -> AddFrame ( _E_minLabel );
  _energy_group_frame -> AddFrame ( _E_min      );
  _energy_group_frame -> AddFrame ( _E_maxLabel );
  _energy_group_frame -> AddFrame ( _E_max      );

  _wcut_group_frame   -> AddFrame ( _W_cut      );

  _cuts_group_frame   -> AddFrame( _energy_group_frame  );
  _cuts_group_frame   -> AddFrame( _wcut_group_frame    );

  //-- bottom/left side: add all parent frames

  _main_left_frame -> AddFrame( _xsec_err_group_frame   );
  _main_left_frame -> AddFrame( _exp_group_frame        );
  _main_left_frame -> AddFrame( _cuts_group_frame       );
  _main_left_frame -> AddFrame( _scale_with_energy      );
}
//______________________________________________________________________________
void vDataSelectionDialog::BuildRightFrameWidgets(void)
{
  _reaction_group_frame   = new TGGroupFrame(_main_right_frame,
                                      "Observable/Reaction",  kHorizontalFrame);
  _obs_group_frame        = new TGGroupFrame(_reaction_group_frame,
                                          "Cross Section",      kVerticalFrame);
  _init_state_group_frame = new TGGroupFrame(_reaction_group_frame,
                                              "Initial State",  kVerticalFrame);
  _target_group_frame     = new TGGroupFrame(_reaction_group_frame,
                                                     "Target",  kVerticalFrame);

  _listbox_layout = new TGLayoutHints(
     kLHintsTop | kLHintsCenterX | kLHintsExpandX | kLHintsExpandY, 1, 1, 1, 1);

  _measurements_listbox = new TGListBox(_main_right_frame,        201);
  _xsec_listbox         = new TGListBox(_obs_group_frame,         202);
  _nu_listbox           = new TGListBox(_init_state_group_frame,  203);
  _target_listbox       = new TGListBox(_target_group_frame,      204);

  _measurements_listbox -> Resize (580, 220);
  _xsec_listbox         -> Resize (100, 90);  
  _nu_listbox           -> Resize (200, 90);
  _target_listbox       -> Resize (170, 90);

  _measurements_listbox -> SetMultipleSelections( true );
  _xsec_listbox         -> SetMultipleSelections( true  );
  _nu_listbox           -> SetMultipleSelections( true  );
  _target_listbox       -> SetMultipleSelections( true  );

  this->LoadMeasurementsFromDB();
  this->LoadXSecTypesFromDB();
  this->LoadProbesFromDB();
  this->LoadTargetsFromDB();

  _select_all_xsec   = new TGCheckButton(_obs_group_frame,        "Select all", 72);
  _select_all_nu     = new TGCheckButton(_init_state_group_frame, "Select all", 73);
  _select_all_target = new TGCheckButton(_target_group_frame,     "Select all", 74);

  _button_layout = new TGLayoutHints( kLHintsTop | kLHintsCenterX, 2, 2, 2, 2 );

  _close_button         = new TGTextButton (_main_right_frame, "&Close", 101);

  
  //--- FORCE CONSISTENCY BETWEEN SELECTED VALUES

  // Update widgets / force consistency every time the "obervables" listbox is clicked

  _xsec_listbox->Connect("SelectionChanged()","genie::nuvld::vDataSelectionDialog",
                                    this, "MakeConsistentWithExpObsListboxes()");

  //--- "Select all" action

  _select_all_xsec   -> Connect("Clicked()","genie::nuvld::vDataSelectionDialog",
                                                        this,"SelectAllXSec()");
  _select_all_nu     -> Connect("Clicked()","genie::nuvld::vDataSelectionDialog",
                                                      this,"SelectAllProbes()");
  _select_all_target -> Connect("Clicked()","genie::nuvld::vDataSelectionDialog",
                                                     this,"SelectAllTargets()");

  _obs_group_frame        -> AddFrame( _xsec_listbox      );
  _obs_group_frame        -> AddFrame( _select_all_xsec   );
  _init_state_group_frame -> AddFrame( _nu_listbox        );
  _init_state_group_frame -> AddFrame( _select_all_nu     );
  _target_group_frame     -> AddFrame( _target_listbox    );
  _target_group_frame     -> AddFrame( _select_all_target );
  
  _reaction_group_frame   -> AddFrame( _obs_group_frame        );
  _reaction_group_frame   -> AddFrame( _init_state_group_frame );
  _reaction_group_frame   -> AddFrame( _target_group_frame     );
                                                    
  _close_button->Connect("Clicked()",
                         "genie::nuvld::vDataSelectionDialog", this, "Close()");  

  _main_right_frame->AddFrame (_measurements_listbox,  _listbox_layout  );
  _main_right_frame->AddFrame (_reaction_group_frame );
  _main_right_frame->AddFrame (_close_button,             _button_layout);
}
//______________________________________________________________________________
void vDataSelectionDialog::LoadExperimentsFromDB(void)
{
  TSQLServer * sql_server = _db->SqlServer();

  const char query[] = "SELECT name from EXP_INFO;";

  TSQLResult * result = sql_server->Query(query);

  const int nrows = result->GetRowCount();

  TSQLRow * row = 0;
    
  for (int i = 0; i < nrows; i++) {

      row = result->Next();

      _exp_listbox->AddEntry( row->GetField(0), i);

      delete row;
  }    
  delete result;

  gClient->NeedRedraw(_exp_listbox->GetContainer());
}
//______________________________________________________________________________
void vDataSelectionDialog::LoadXSecTypesFromDB(void)
{
  TSQLServer * sql_server = _db->SqlServer();

  const char query[] =
      "SELECT observable from MEASUREMENT_HEADER WHERE reaction LIKE \"%nu%\";";

  TSQLResult * result = sql_server->Query(query);

  map<string, int> xsec_types;

  const int nrows = result->GetRowCount();

  TSQLRow * row = 0;

  for (int i = 0; i < nrows; i++) {

      row = result->Next();

      // add in a map to remove duplicate keys      
      xsec_types[string(row->GetField(0))]++; 

      delete row;
  }
  delete result;

  map<string, int>::const_iterator xsec_types_iter;

  int i=0;
  
  for(xsec_types_iter = xsec_types.begin();
         xsec_types_iter != xsec_types.end(); ++xsec_types_iter)
                  _xsec_listbox->AddEntry( xsec_types_iter->first.c_str(), i++);

  gClient->NeedRedraw(_xsec_listbox->GetContainer());                  
}
//______________________________________________________________________________
void vDataSelectionDialog::LoadProbesFromDB(void)
{
  TSQLServer * sql_server = _db->SqlServer();

  const char query[] =
      "SELECT reaction from MEASUREMENT_HEADER WHERE reaction LIKE \"%nu%\";";

  TSQLResult * result = sql_server->Query(query);

  map<string, int> xsec_types;

  const int nrows = result->GetRowCount();

  TSQLRow * row = 0;

  for (int i = 0; i < nrows; i++) {

      row = result->Next();

      // add in a map to remove duplicate keys
      xsec_types[string(row->GetField(0))]++;

      delete row;
  }
  delete result;

  map<string, int>::const_iterator xsec_types_iter;

  int i=0;

  for(xsec_types_iter = xsec_types.begin();
         xsec_types_iter != xsec_types.end(); ++xsec_types_iter)
                  _nu_listbox->AddEntry( xsec_types_iter->first.c_str(), i++);

  gClient->NeedRedraw(_nu_listbox->GetContainer());                  
}
//______________________________________________________________________________
void vDataSelectionDialog::LoadTargetsFromDB(void)
{
  TSQLServer * sql_server = _db->SqlServer();

  const char query[] =
         "SELECT target from MEASUREMENT_HEADER WHERE reaction LIKE \"%nu%\";";

  TSQLResult * result = sql_server->Query(query);

  map<string, int> targets;

  const int nrows = result->GetRowCount();

  TSQLRow * row = 0;

  for (int i = 0; i < nrows; i++) {

      row = result->Next();

      // add in a map to remove duplicate keys
      targets[string(row->GetField(0))]++;

      delete row;
  }
  delete result;

  map<string, int>::const_iterator target_iter;

  int i=0;

  for(target_iter = targets.begin();
             target_iter != targets.end(); ++target_iter)
                    _target_listbox->AddEntry( target_iter->first.c_str(), i++);

  gClient->NeedRedraw(_target_listbox->GetContainer());                    
}
//______________________________________________________________________________
void vDataSelectionDialog::LoadMeasurementsFromDB(void)
{
  TSQLServer * sql_server = _db->SqlServer();

  const char query[] =
        "SELECT MEASUREMENT_HEADER.observable, MEASUREMENT_HEADER.reaction, \
         MEASUREMENT_HEADER.name, MEASUREMENT_HEADER.measurement_tag, \
         REFERENCE.authors, REFERENCE.journal, REFERENCE.year \
         FROM MEASUREMENT_HEADER, REFERENCE \
         WHERE REFERENCE.name = MEASUREMENT_HEADER.name AND \
         REFERENCE.measurement_tag = MEASUREMENT_HEADER.measurement_tag \
         AND MEASUREMENT_HEADER.reaction LIKE \"%nu%\";";

  TSQLResult * result = sql_server->Query(query);

  const int nrows = result->GetRowCount();

  TSQLRow * row = 0;

  for (int i = 0; i < nrows; i++) {

      row = result->Next();

      ostringstream item;

      item << row->GetField(2) << ";"
           << row->GetField(3) << ";"
           << row->GetField(4) << ";"
           << row->GetField(5) << ";"
           << row->GetField(6) << ";"
           << row->GetField(0) << ";"
           << row->GetField(1);
                      
      _measurements_listbox->AddEntry(item.str().c_str(), i);

      cout << item.str() << endl;

      delete row;
  }
  delete result;

  gClient->NeedRedraw(_measurements_listbox->GetContainer());
}
//______________________________________________________________________________
void vDataSelectionDialog::MakeConsistentWithExpListbox(void)
{
  TList * selected = new TList();
  
  _exp_listbox->GetSelectedEntries(selected);

  TGTextLBEntry * selected_entry = 0;

  TIter selected_iter(selected);

  //--- reset all dependent:
  
  gui_utils::ResetAllListBoxSelections(_measurements_listbox); // measurement list entries 
  gui_utils::ResetAllListBoxSelections(_xsec_listbox); // observables
  gui_utils::ResetAllListBoxSelections(_nu_listbox); // reactions
  gui_utils::ResetAllListBoxSelections(_target_listbox); // targets
  
  while( (selected_entry = (TGTextLBEntry *) selected_iter.Next()) ) {
  
     string exp_name = selected_entry->GetText()->GetString();

     //--- select al matched "measurement list" entries
     SelectAllMeasurementsForExp(exp_name);

     //--- select all matched "observables"      
     SelectAllObservablesForExp(exp_name);

     //--- select all matched "reactions" 
     SelectAllReactionsForExp(exp_name);

     //--- select all matched "targets" 
     SelectAllTargetsForExp(exp_name);

  }
  delete selected;
}
//______________________________________________________________________________
void vDataSelectionDialog::MakeConsistentWithExpObsListboxes(void)
{
  TList * selected_exp = new TList();
  TList * selected_obs = new TList();

  _exp_listbox  -> GetSelectedEntries(selected_exp);
  _xsec_listbox -> GetSelectedEntries(selected_obs);

  TGTextLBEntry * exp_entry = 0;
  TGTextLBEntry * obs_entry = 0;

  TIter exp_iter(selected_exp);
  TIter obs_iter(selected_obs);

  //--- reset all dependent:

  gui_utils::ResetAllListBoxSelections(_measurements_listbox);
  gui_utils::ResetAllListBoxSelections(_nu_listbox); 
  gui_utils::ResetAllListBoxSelections(_target_listbox); 

  while( (exp_entry = (TGTextLBEntry *) exp_iter.Next()) ) {

     string exp_name = exp_entry->GetText()->GetString();

     obs_iter.Reset();
     
     while( (obs_entry = (TGTextLBEntry *) obs_iter.Next()) ) {

         string observable = obs_entry->GetText()->GetString();

         //--- select al matched "measurement list" entries
         SelectAllMeasurementsForExpObs(exp_name, observable);

         //--- select all matched "reactions"
         SelectAllReactionsForExpObs(exp_name, observable);

         //--- select all matched "targets"
         SelectAllTargetsForExpObs(exp_name, observable);
     }         
  }
  
  delete selected_exp;
  delete selected_obs;
}
//______________________________________________________________________________
void vDataSelectionDialog::SelectAllMeasurementsForExp(string exp_name)
{
  //--- measurement list entries matching with a selected experiment name
  //    are set true

  int nmle = _measurements_listbox->GetNumberOfEntries();

  for(int i = 0; i < nmle; i++) {

     TGTextLBEntry* mle = (TGTextLBEntry*) _measurements_listbox->GetEntry(i);

     string mle_as_string = mle->GetText()->GetString();

     if(mle_as_string.find(exp_name) != string::npos)
                                        _measurements_listbox->Select(i, true);                                        
  }
  
  gClient->NeedRedraw(_measurements_listbox->GetContainer());
}
//______________________________________________________________________________
void vDataSelectionDialog::SelectAllObservablesForExp(string exp_name)
{
  TSQLRow * row = 0;

  ostringstream query;

  query << "SELECT observable FROM MEASUREMENT_HEADER where name = \""
        << exp_name << "\";";

  TSQLResult * result = _db->SqlServer()->Query( query.str().c_str() );

  const int nobs  = _xsec_listbox->GetNumberOfEntries();

  const int nrows = result->GetRowCount();

  for (int i = 0; i < nrows; i++) {

     row = result->Next();

     string exp_obs = row->GetField(0); // current observable

     // loop over the observables listbox and highlight the matched entries
     
     for(int i = 0; i < nobs; i++) {

        TGTextLBEntry* obs = (TGTextLBEntry*) _xsec_listbox->GetEntry(i);

        string obs_as_string = obs->GetText()->GetString();

        if(obs_as_string.find(exp_obs) != string::npos)
                                                _xsec_listbox->Select(i, true);
     }                                              
  }
  delete result;

  gClient->NeedRedraw(_xsec_listbox->GetContainer());  
}
//______________________________________________________________________________
void vDataSelectionDialog::SelectAllReactionsForExp(string exp_name)
{
  TSQLRow * row = 0;

  ostringstream query;

  query << "SELECT reaction FROM MEASUREMENT_HEADER where name = \""
        << exp_name << "\";";

  TSQLResult * result = _db->SqlServer()->Query( query.str().c_str() );

  const int nreac  = _nu_listbox->GetNumberOfEntries();

  const int nrows = result->GetRowCount();

  for (int i = 0; i < nrows; i++) {

     row = result->Next();

     string exp_reac = row->GetField(0); // current reaction

     // loop over the reactions listbox and highlight the matched entries

     for(int i = 0; i < nreac; i++) {

        TGTextLBEntry* reac = (TGTextLBEntry*) _nu_listbox->GetEntry(i);

        string reac_as_string = reac->GetText()->GetString();

        if(reac_as_string.find(exp_reac) != string::npos)
                                                _nu_listbox->Select(i, true);
     }
  }
  delete result;

  gClient->NeedRedraw(_nu_listbox->GetContainer());  
}
//______________________________________________________________________________
void vDataSelectionDialog::SelectAllTargetsForExp(string exp_name)
{
  TSQLRow * row = 0;

  ostringstream query;

  query << "SELECT target FROM MEASUREMENT_HEADER where name = \""
        << exp_name << "\";";

  TSQLResult * result = _db->SqlServer()->Query( query.str().c_str() );

  const int ntargets  = _target_listbox->GetNumberOfEntries();

  const int nrows = result->GetRowCount();

  for (int i = 0; i < nrows; i++) {

     row = result->Next();

     string exp_target = row->GetField(0); // current target

     // loop over the targets listbox and highlight the matched entries

     for(int i = 0; i < ntargets; i++) {

        TGTextLBEntry* target = (TGTextLBEntry*) _target_listbox->GetEntry(i);

        string target_as_string = target->GetText()->GetString();

        if(target_as_string.find(exp_target) != string::npos)
                                                _target_listbox->Select(i, true);
     }
  }
  delete result;

  gClient->NeedRedraw(_target_listbox->GetContainer());  
}
//______________________________________________________________________________
void vDataSelectionDialog::SelectAllMeasurementsForExpObs(
                                             string exp_name, string observable)
{
  //--- measurement list entries matching with a selected experiment name
  //    and observable are set true

  int nmle = _measurements_listbox->GetNumberOfEntries();

  for(int i = 0; i < nmle; i++) {

     TGTextLBEntry* mle = (TGTextLBEntry*) _measurements_listbox->GetEntry(i);

     string mle_as_string = mle->GetText()->GetString();

     if(mle_as_string.find(exp_name) != string::npos &&
                   mle_as_string.find(observable) != string::npos)
                                        _measurements_listbox->Select(i, true);
  }

  gClient->NeedRedraw(_measurements_listbox->GetContainer());
}
//______________________________________________________________________________
void vDataSelectionDialog::SelectAllReactionsForExpObs(
                                            string exp_name, string observable)
{
  TSQLRow * row = 0;

  ostringstream query;

  query << "SELECT reaction FROM MEASUREMENT_HEADER where "
        << "name = \"" << exp_name << "\" AND "
        << "observable = \"" << observable << "\";";

  TSQLResult * result = _db->SqlServer()->Query( query.str().c_str() );

  const int nreac  = _nu_listbox->GetNumberOfEntries();

  const int nrows = result->GetRowCount();

  for (int i = 0; i < nrows; i++) {

     row = result->Next();

     string exp_reac = row->GetField(0); // current reaction

     // loop over the reactions listbox and highlight the matched entries

     for(int i = 0; i < nreac; i++) {

        TGTextLBEntry* reac = (TGTextLBEntry*) _nu_listbox->GetEntry(i);

        string reac_as_string = reac->GetText()->GetString();

        if(reac_as_string.find(exp_reac) != string::npos)
                                                _nu_listbox->Select(i, true);
     }
  }
  delete result;

  gClient->NeedRedraw(_nu_listbox->GetContainer());
}
//______________________________________________________________________________
void vDataSelectionDialog::SelectAllTargetsForExpObs(
                                             string exp_name, string observable)
{
  TSQLRow * row = 0;

  ostringstream query;

  query << "SELECT target FROM MEASUREMENT_HEADER where "
        << "name = \"" << exp_name << "\" AND "
        << "observable = \"" << observable << "\";";

  TSQLResult * result = _db->SqlServer()->Query( query.str().c_str() );

  const int ntargets  = _target_listbox->GetNumberOfEntries();

  const int nrows = result->GetRowCount();

  for (int i = 0; i < nrows; i++) {

     row = result->Next();

     string exp_target = row->GetField(0); // current target

     // loop over the targets listbox and highlight the matched entries

     for(int i = 0; i < ntargets; i++) {

        TGTextLBEntry* target = (TGTextLBEntry*) _target_listbox->GetEntry(i);

        string target_as_string = target->GetText()->GetString();

        if(target_as_string.find(exp_target) != string::npos)
                                                _target_listbox->Select(i, true);
     }
  }
  delete result;

  gClient->NeedRedraw(_target_listbox->GetContainer());
}
//______________________________________________________________________________
void vDataSelectionDialog::SelectAllExp(void)
{
  if(_select_all_exp->GetState() == kButtonDown)
                               gui_utils::SelectAllListBoxEntries(_exp_listbox);
  else gui_utils::ResetAllListBoxSelections(_exp_listbox);

  _exp_listbox->Resize(120, 90);

  gClient->NeedRedraw(_exp_listbox->GetContainer());

  MakeConsistentWithExpListbox();
}
//______________________________________________________________________________
void vDataSelectionDialog::SelectAllXSec(void)
{
  if(_select_all_xsec->GetState() == kButtonDown)
                              gui_utils::SelectAllListBoxEntries(_xsec_listbox);
  else gui_utils::ResetAllListBoxSelections(_xsec_listbox);

  gClient->NeedRedraw(_xsec_listbox->GetContainer());  
}
//______________________________________________________________________________
void vDataSelectionDialog::SelectAllProbes(void)
{
  if(_select_all_nu->GetState() == kButtonDown)
                                gui_utils::SelectAllListBoxEntries(_nu_listbox);
  else gui_utils::ResetAllListBoxSelections(_nu_listbox);

  gClient->NeedRedraw(_nu_listbox->GetContainer());    
}
//______________________________________________________________________________
void vDataSelectionDialog::SelectAllTargets(void)
{
  if(_select_all_target->GetState() == kButtonDown)
                            gui_utils::SelectAllListBoxEntries(_target_listbox);
  else gui_utils::ResetAllListBoxSelections(_target_listbox);

  gClient->NeedRedraw(_target_listbox->GetContainer());  
}
//______________________________________________________________________________
string vDataSelectionDialog::BundleKeyListInString(void)
{
  unsigned int ikey = 0;

  ostringstream key_list;

  TList * selected = new TList();

  _measurements_listbox->GetSelectedEntries(selected);

  TGTextLBEntry * entry = 0;

  TIter iter(selected);

  int nselected = selected->IndexOf( selected->Last() ) + 1;
  
  while( (entry = (TGTextLBEntry *) iter.Next()) ) {

      vector<string> key_elem = ParserUtils::split(
                                     entry->GetText()->GetString(),  ";");

      assert( key_elem.size() == 7 );

      key_list << key_elem[0] << "," << key_elem[1];

      if(ikey++ < nselected-1) key_list << ";";
  }

  return key_list.str().c_str();
}
//______________________________________________________________________________
string vDataSelectionDialog::BundleCutsInString(void)
{
  float Emin = _E_min->GetNumber();
  float Emax = _E_max->GetNumber();
  
  ostringstream cuts;

  cuts << "Emin=" << Emin << ";" << "Emax=" << Emax;

  return cuts.str();
}
//______________________________________________________________________________
string vDataSelectionDialog::BundleDrawOptInString(void)
{
  if(_scale_with_energy->GetState() == kButtonDown) return "scale-with-energy";
  else return "";
}
//______________________________________________________________________________  
void vDataSelectionDialog::PositionRelativeToParent(const TGWindow * main)
{
// position relative to the parent's window
  
  int ax, ay;
  Window_t wdum;

  gVirtualX->TranslateCoordinates(
              main->GetId(), _main->GetParent()->GetId(),
              (Int_t)(((TGFrame *) main)->GetWidth() - _main->GetWidth()) >> 1,
              (Int_t)(((TGFrame *) main)->GetHeight() - _main->GetHeight()) >> 1,
              ax, ay, wdum);
              
  _main->Move(ax, ay);
}
//______________________________________________________________________________


