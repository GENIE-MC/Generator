//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiStackHandler

\brief    Responds to GUI events associated with stacking and retrieving from
          the stack of DBTable<T>'s and generator configurations.

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 26, 2004
*/
//_____________________________________________________________________________

#ifndef _GUI_STACK_HANDLER_H_
#define _GUI_STACK_HANDLER_H_

#include <vector>
#include <string>

#include <TGComboBox.h>
#include <TGTextEntry.h>

#include  "ValidationTools/NuVld/DBI.h"
#include  "ValidationTools/NuVld/DBTable.h"
#include  "ValidationTools/NuVld/DBTableStack.h"
#include  "ValidationTools/NuVld/DBNuXSecTableRow.h"
#include  "ValidationTools/NuVld/DBElDiffXSecTableRow.h"
#include  "ValidationTools/NuVld/DBConnection.h"

using std::vector;
using std::string;

namespace genie {
namespace nuvld {

class GuiStackHandler {

public:
  GuiStackHandler();
 ~GuiStackHandler();

  //-- attach the 'data/config stacking' GUI widgets
  
  void SetParentMainFrame    (TGMainFrame * mf) { fMain = mf; }

  void AttachDBTableTxtEntry (TGTextEntry * te) { fDBTableTxtEntry     = te; }
  void AttachConfigTxtEntry  (TGTextEntry * te) { fConfigTxtEntry      = te; }
  void AttachDBTableCombo    (TGComboBox *  cb) { fStackedDBTableCombo = cb; }
  void AttachConfigCombo     (TGComboBox *  cb) { fStackedConfigCombo  = cb; }

  void SetDBConnection (DBConnection * dbc) { fDBC = dbc; }
  
  //-- response to 'stack GUI' events
  
  void SaveStack                 (void);
  void LoadStack                 (void);
  void StackDBTable              (void);
  void StackConfig               (void);  
  void EraseStackedItem          (void);

  //-- methods for saving / loading the stack to / from ROOT files
  
  void LoadDBTables       (string root_file, bool keep_current_stack);
  void LoadNeugenConfig   (string root_file, bool keep_current_stack);

  //-- methods for updating the combo-boxes once some entry is added/deleted

  void   UpdateStackedDBTableCombo (void);
  void   UpdateStackedConfigCombo  (void); 
  void   EraseStackedDBTable       (void);
  void   EraseStackedConfig        (void);
  string StackedDBTableName        (unsigned int id) const;
  string StackedConfigName         (unsigned int id) const;
    
private:

  TGMainFrame *  fMain;

  TGTextEntry *  fDBTableTxtEntry; 
  TGTextEntry *  fConfigTxtEntry;
  TGComboBox *   fStackedDBTableCombo;
  TGComboBox *   fStackedConfigCombo;

  DBConnection * fDBC;
};

} // nuvld namespace
} // genie namespace

#endif

