//_____________________________________________________________________________
/*!

\class    genie::nuvld::DataSelectionDialog

\brief    Base class for data selection popup dialogs and tabs

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  November 01, 2004
*/
//_____________________________________________________________________________

#include <sstream>
#include <iostream>

#include "NuVldGUI/DataSelectionDialog.h"

using std::ostringstream;

using namespace genie::nuvld;

//______________________________________________________________________________
DataSelectionDialog::DataSelectionDialog(void) 
{

}
//______________________________________________________________________________
DataSelectionDialog::~DataSelectionDialog()
{

}
//______________________________________________________________________________
string DataSelectionDialog::BundleSelectionsInString(void)
{
  ostringstream options;

  options << "KEY-LIST:" << BundleKeyListInString() << "$"
          << "CUTS:"     << BundleCutsInString()    << "$"
          << "DRAW_OPT:" << BundleDrawOptInString() << "$"
          << "DB-TYPE:vN-XSec"; // tmp

  return options.str();
}
//______________________________________________________________________________

