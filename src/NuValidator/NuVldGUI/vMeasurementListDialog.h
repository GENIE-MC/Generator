//_____________________________________________________________________________
/*!

\class    genie::nuvld::vMeasurementListDialog

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  November 01, 2004
*/
//_____________________________________________________________________________

#ifndef _NEUTRINO_MEASUREMENT_LIST_DIALOG_H_
#define _NEUTRINO_MEASUREMENT_LIST_DIALOG_H_

#include <TApplication.h>
#include <TVirtualX.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGListBox.h>
#include <TGButton.h>
#include <RQ_OBJECT.h>

#include "NuVldGUI/DBConnection.h"
#include "NuVldGUI/DataSelectionDialog.h"

namespace genie {
namespace nuvld {

class vMeasurementListDialog : public DataSelectionDialog {

RQ_OBJECT("vMeasurementListDialog")

public:
   vMeasurementListDialog(
            const TGWindow *p, const TGWindow *main, bool & attn,
            UInt_t w, UInt_t h, UInt_t options = kVerticalFrame, DBConnection * db = 0);
   virtual ~vMeasurementListDialog();

   void CloseWindow (void) { delete this;               }
   void Close       (void) { _main->SendCloseMessage(); }

   //--- DataSelectionDialog interface
   
   string BundleKeyListInString (void);
   string BundleCutsInString    (void);
   string BundleDrawOptInString (void);
   
private:

   void LoadMeasurementsFromDB   (void);
   void PositionRelativeToParent (const TGWindow * main);

   TGTransientFrame * _main;
   TGListBox *        _measurements_listbox;
   TGButton *         _close_button;
   TGLayoutHints *    _listbox_layout;
   TGLayoutHints *    _button_layout;

   DBConnection *     _db;

   ClassDef(vMeasurementListDialog, 0)
};

} // nuvld namespace
} // genie namespace

#endif // _NEUTRINO_MEASUREMENT_LIST_DIALOG_H_

