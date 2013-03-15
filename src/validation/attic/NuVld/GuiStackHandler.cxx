//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Aug 26, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include <sstream>

#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TGFileDialog.h>
#include <TGProgressBar.h>

#include "Messenger/Messenger.h"
#include "Utils/GUIUtils.h"
#include "ValidationTools/NuVld/GuiStackHandler.h"
#include "ValidationTools/NuVld/NuVldUserData.h"
#include "ValidationTools/NuVld/GuiSysLogSingleton.h"
#include "ValidationTools/NuVld/GuiYNQuestionBox.h"
#include "ValidationTools/NuVld/GuiMsgBox.h"

using std::ostringstream;

using namespace genie;
using namespace genie::nuvld;

//______________________________________________________________________________
GuiStackHandler::GuiStackHandler()
{
  fMain = 0;

  fDBTableTxtEntry     = 0;
  fConfigTxtEntry      = 0;
  fStackedDBTableCombo = 0;
//fStackedConfigCombo  = 0;

  fDBC = 0;
}
//______________________________________________________________________________
GuiStackHandler::~GuiStackHandler()
{

}
//______________________________________________________________________________
void GuiStackHandler::SaveStack(void)
{
  GuiSysLogSingleton * syslog  = GuiSysLogSingleton::Instance();
  NuVldUserData * user_data = NuVldUserData::Instance();

  static TString dir(".");

  const char * kFileExt[] = { "All files", "*", "ROOT files", "*.root", 0, 0 };

  TGFileInfo fi;
  fi.fFileTypes = kFileExt;
  fi.fIniDir    = StrDup( dir.Data() );

  new TGFileDialog(gClient->GetRoot(), fMain, kFDSave, &fi);

  if( fi.fFilename ) {

     string root_file = string( fi.fFilename );

     ostringstream cmd;
     cmd << "Saving stacked data to file: " << root_file.c_str();

     syslog -> Log()       -> AddLine( cmd.str().c_str()    );
     syslog -> StatusBar() -> SetText( cmd.str().c_str(), 0 );

     user_data->SaveStack(root_file);

     // fancy fake effect...
     for(int i=0; i<100; i++) {
        gSystem->Sleep(3);
        syslog -> ProgressBar() ->SetPosition(i);
     }
     syslog -> ProgressBar() ->SetPosition(0);
  }
}
//______________________________________________________________________________
void GuiStackHandler::LoadStack(void)
{
  NuVldUserData * user_data = NuVldUserData::Instance();
  GuiSysLogSingleton * syslog  = GuiSysLogSingleton::Instance();

  // check if there are already stacked datasets & configs from this session
  // and then ask if he wants to keep them together with the items he loads

  bool keep_current_stacked_data = true;

  const char * msg_data =
              " You already have some stacked data-sets. Should I keep them? ";

  if( user_data->NStackedDBTables() > 0 ) {

     new GuiYNQuestionBox(gClient->GetRoot(),
         fMain, 380, 250, kVerticalFrame, msg_data, &keep_current_stacked_data);
  }

/*
  bool keep_current_stacked_configs = true;

  const char * msg_conf =
         " You already have some stacked NeuGEN configs. Should I keep them? ";

  if( user_data->NStackedConfigs() > 0 ) {

       new GuiYNQuestionBox(gClient->GetRoot(),
      fMain, 380, 250, kVerticalFrame, msg_conf, &keep_current_stacked_configs);
  }
*/
  // root file-open dialog

  static TString dir(".");

  const char * kFileExt[] = { "All files", "*", "ROOT files", "*.root", 0, 0 };

  TGFileInfo fi;
  fi.fFileTypes = kFileExt;
  fi.fIniDir    = StrDup( dir.Data() );

  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);

  if( fi.fFilename ) {

     // get selected filename

     string root_file = string( fi.fFilename );

     // session log & status bar messages

     ostringstream cmd;
     cmd << "Loading stacked data from file: " << root_file.c_str();

     syslog -> Log()       -> AddLine( cmd.str().c_str()    );
     syslog -> StatusBar() -> SetText( cmd.str().c_str(), 0 );

   // this -> LoadNeugenConfig (root_file, keep_current_stacked_configs);
     this -> LoadDBTables     (root_file, keep_current_stacked_data);

     // update the combo box with all the stacked data items

     this -> UpdateStackedDBTableCombo ();
   //this -> UpdateStackedConfigCombo  ();

     // fancy fake effect...
     for(int i=0; i<100; i++) {
        gSystem->Sleep(3);
        syslog -> ProgressBar() ->SetPosition(i);
     }
     syslog -> ProgressBar() -> SetPosition(0);
  }
}
//______________________________________________________________________________
void GuiStackHandler::StackDBTable(void)
{
  LOG("NuVld", pDEBUG) << "Attempting to stack a DBTable<T>";

  NuVldUserData * user_data = NuVldUserData::Instance();

  // get the contents of the TGTextEntry

  string name = fDBTableTxtEntry->GetBuffer()->GetString();

  LOG("NuVld", pDEBUG) << "Name given by the user to the DBTable<T>: " << name;

  // check whether the given name or current DBTable<T> is null

  if( user_data->CurrDBTableIsNull() )
      new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                                       "May I notice that YOU HAVE NO DATA?");

  else if ( ! (name.size() > 0) )
      new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                       "Didn't you forget to give a name to your selection?");
  else {

     LOG("NuVld", pDEBUG)
                    << "Asking NuVldUserData to stack the current DBTable<T>";

     user_data->AddCurrDBTableToStack(name);

     // update the combo box with all the stacked data items

     this->UpdateStackedDBTableCombo();

     // reset the TGTextEntry

     fDBTableTxtEntry->GetBuffer()->Clear();
     fDBTableTxtEntry->SetCursorPosition(0);

     gClient->NeedRedraw(fDBTableTxtEntry); // force re-fresh
  }
}
//______________________________________________________________________________
void GuiStackHandler::StackConfig(void)
{
  LOG("NuVld", pDEBUG) << "Attempting to stack a generator configuration";
/*
  NuVldUserData * user_data = NuVldUserData::Instance();

  // get the contents of the TGTextEntry

  string name = fConfigTxtEntry->GetBuffer()->GetString();

  LOG("NuVld", pDEBUG) << "Name given by the user to the config.: " << name;

  // if there is something written there

  if( name.size() > 0 ) {

     // add the new named pair of cards

     NeuGenCards * cards = NeuGenCards::Instance();

     NGCardPair * pair = new NGCardPair;

     pair->SetConfig( cards->CurrConfig() );
     pair->SetInputs( cards->CurrInputs() );

     LOG("NuVld", pDEBUG)
                    << "Asking NuVldUserData to stack the current config";

     user_data -> NeuGenCardPairStack() -> AddCardPair(name, pair);

     // update the corresponding TComboBox with the stacked entries

     this -> UpdateStackedConfigCombo ();

     // reset the TGTextEntry

     fConfigTxtEntry->GetBuffer()->Clear();
     fConfigTxtEntry->SetCursorPosition(0);

     gClient->NeedRedraw(fConfigTxtEntry); // force re-fresh

    } else
       new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                                               "You must enter a name first!");
*/
}
//______________________________________________________________________________
void GuiStackHandler::EraseStackedItem(void)
{
  if(fStackedDBTableCombo->GetSelectedEntry()) this->EraseStackedDBTable();
//if(fStackedConfigCombo->GetSelectedEntry() ) this->EraseStackedConfig ();
}
//______________________________________________________________________________
void GuiStackHandler::EraseStackedDBTable(void)
{
  NuVldUserData * user_data = NuVldUserData::Instance();

  if(fStackedDBTableCombo->GetSelectedEntry()) {

    int id = fStackedDBTableCombo->GetSelectedEntry()->EntryId();

    if(id >= 0 && id < (int) user_data->NStackedDBTables()) {

       string dbtname = this->StackedDBTableName((unsigned int)id);

       user_data->DelStackedDBTable(dbtname);

       // update the combo box with all the stacked data items

       this->UpdateStackedDBTableCombo();

    }//entry-id

  }//selected-entry
}
//______________________________________________________________________________
void GuiStackHandler::EraseStackedConfig(void)
{
/*
  NuVldUserData * user_data = NuVldUserData::Instance();

  if(fStackedConfigCombo->GetSelectedEntry()) {

    int id = fStackedConfigCombo->GetSelectedEntry()->EntryId();

    if(id >= 0 && id < (int) user_data->NStackedConfigs()) {

       // delete the selected config

       string cfgname = this->StackedConfigName((unsigned int)id);

       user_data->NeuGenCardPairStack()->Erase(cfgname);

       // update the corresponding TComboBox with the stacked entries

       this -> UpdateStackedConfigCombo ();

    }//entry-id
  }//selected-entry
*/
}
//______________________________________________________________________________
void GuiStackHandler::LoadDBTables(string root_file, bool keep_current)
{
  NuVldUserData * user_data = NuVldUserData::Instance();

  if( ! keep_current)
  {
     LOG("NuVld", pDEBUG)
            << "Asking NuVldUserData to clear existing stacked DBTable<T>'s";

     user_data->ClearStackedDBTables();
  }

  bool is_connected = (fDBC->SqlServer() != 0 && fDBC->SqlServer()->IsConnected());

  if( is_connected ) {

    LOG("NuVld", pDEBUG) << "Found an active DB connection";

    DBI dbi( fDBC->SqlServer() );

    TFile f(root_file.c_str(),"READ");

    TDirectory * v_xsec_dir = (TDirectory *) f.Get("v_xsec");

    if(v_xsec_dir) {

       TKey *  key = 0;
       TList * keys = v_xsec_dir->GetListOfKeys();
       TIter   key_iter(keys);

       while ( (key = (TKey *) key_iter.Next()) ) {

           DBTable<DBNuXSecTableRow> * table = new DBTable<DBNuXSecTableRow>;

           LOG("NuVld", pDEBUG) << "Loading:..........." << key->GetName();

           DBQueryString * query_string =
                          (DBQueryString *) v_xsec_dir->Get(key->GetName());

           dbi.FillTable(table, *query_string);

           user_data->NuXSecStack()->AddDBTable( key->GetName(), table);
       }
    }

    TDirectory * e_diff_xsec_dir = (TDirectory *) f.Get("e_diff_xsec");

    if(e_diff_xsec_dir) {

       TKey *  key  = 0;
       TList * keys = e_diff_xsec_dir->GetListOfKeys();
       TIter   key_iter(keys);

       while ( (key = (TKey *) key_iter.Next()) ) {

           DBTable<DBElDiffXSecTableRow> * table = new DBTable<DBElDiffXSecTableRow>;

           LOG("NuVld", pDEBUG) << "Loading:..........." << key->GetName();

           DBQueryString * query_string =
                         (DBQueryString *) e_diff_xsec_dir->Get(key->GetName());

           dbi.FillTable(table, *query_string);

           user_data->ElDiffXSecStack()->AddDBTable( key->GetName(), table);
       }
    }

    f.Close();

  } else {

     LOG("NuVld", pERROR)
              << "Can not load saved stack without connection to the RDBMS";

     new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                                    "You must connect to the data-base first!");

  }
}
//______________________________________________________________________________
void GuiStackHandler::LoadNeugenConfig(
   string /*root_file*/, bool /*keep_current*/)
{
/*
  TFile f(root_file.c_str(), "READ");

  NGCardPairList * cp = (NGCardPairList *) f.Get("generator_config");

  NuVldUserData * user_data = NuVldUserData::Instance();

  // either make this a new _card_pairs or append it to the existing one

  if( ! keep_current)
  {
     LOG("NuVld", pDEBUG)
            << "Asking NuVldUserData to clear existing stacked configurations";

     user_data->ClearStackedConfig();
  }

  if(cp)
  {
     LOG("NuVld", pDEBUG) << "Adding loaded configurations";

     user_data->NeuGenCardPairStack()->Merge(cp);
  }
  f.Close();
*/
}
//______________________________________________________________________________
void GuiStackHandler::UpdateStackedDBTableCombo(void)
{
  LOG("NuVld", pDEBUG) << "Updating Stacked DBTable<T> combo-box";

  NuVldUserData * user_data = NuVldUserData::Instance();

  int nentries = fStackedDBTableCombo->GetNumberOfEntries();

  LOG("NuVld", pDEBUG) << "Removing combo-box entries: 0 --> " << nentries;

  fStackedDBTableCombo->RemoveEntries(0,nentries);

  const vector<string> * names = user_data->GetStackedDBTableNames();

  if(names) {
    utils::gui::FillComboBox( fStackedDBTableCombo, names );

    fStackedDBTableCombo->Select(0);
    gClient->ForceRedraw();

    delete names;
  }
}
//______________________________________________________________________________
void GuiStackHandler::UpdateStackedConfigCombo(void)
{
/*
  LOG("NuVld", pDEBUG) << "Updating Stacked Configurations combo-box";

  NuVldUserData * user_data = NuVldUserData::Instance();

  int nentries = fStackedConfigCombo->GetNumberOfEntries();

  LOG("NuVld", pDEBUG) << "Removing combo-box entries: 0 --> " << nentries;

  fStackedConfigCombo->RemoveEntries(0,nentries);

  const vector<string> * names =
                           user_data->NeuGenCardPairStack()->GetListOfNames();

  if(names) {
    utils::gui::FillComboBox( fStackedConfigCombo, names );

    fStackedConfigCombo->Select(0);
    gClient->ForceRedraw();

    delete names;
  }
*/
}
//______________________________________________________________________________
string GuiStackHandler::StackedDBTableName(unsigned int id) const
{
  NuVldUserData * user_data = NuVldUserData::Instance();

  const vector<string> * names = user_data->GetStackedDBTableNames();

  string dbtname = "";

  if(names) {
     if(id<names->size()) dbtname = (*names)[id];
     delete names;
  }

  return dbtname;
}
//______________________________________________________________________________
string GuiStackHandler::StackedConfigName(unsigned int /*id*/) const
{
/*
  NuVldUserData * user_data = NuVldUserData::Instance();

  const vector<string> * names =
                            user_data->NeuGenCardPairStack()->GetListOfNames();
  string cfgname = "";

  if(names) {
     if(id<names->size()) cfgname = (*names)[id];
     delete names;
  }

  return cfgname;
*/
  return "";
}
//______________________________________________________________________________
