//_____________________________________________________________________________
/*!

\class    genie::nuvld::DataSelectionDialog

\brief    Base class for data selection popup dialogs and tabs

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  November 01, 2004
*/
//_____________________________________________________________________________

#ifndef _DATA_SELECTION_DIALOG_H_
#define _DATA_SELECTION_DIALOG_H_

#include <string>

using std::string;

namespace genie {
namespace nuvld {

class DataSelectionDialog {

public:

   //--- data selection dialog interface
   
   virtual string BundleSelectionsInString (void);   
   virtual string BundleKeyListInString    (void) = 0;
   virtual string BundleCutsInString       (void) = 0;
   virtual string BundleDrawOptInString    (void) = 0;
   virtual void   ResetSelections          (void) = 0;

protected:

   DataSelectionDialog();
   virtual ~DataSelectionDialog();
};

} // nuvld namespace
} // genie namespace

#endif // _DATA_SELECTION_DIALOG_H_

