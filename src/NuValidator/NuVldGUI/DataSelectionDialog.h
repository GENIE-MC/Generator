//_____________________________________________________________________________
/*!

\class    genie::nuvld::DataSelectionDialog

\brief    Base class for data selection pop-up dialogs that require (& lock) 
          the attention of the main GUI throughout their lifetime

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
   
   string  BundleSelectionsInString(void);
   
   virtual string BundleKeyListInString (void) = 0;
   virtual string BundleCutsInString    (void) = 0;
   virtual string BundleDrawOptInString (void) = 0;

protected:

   DataSelectionDialog(bool & attn);
   virtual ~DataSelectionDialog();

   bool & _attn;
};

} // nuvld namespace
} // genie namespace

#endif // _DATA_SELECTION_DIALOG_H_

