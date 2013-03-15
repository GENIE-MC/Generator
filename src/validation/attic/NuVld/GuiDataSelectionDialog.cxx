//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include <sstream>
#include <iostream>

#include "ValidationTools/NuVld/GuiDataSelectionDialog.h"

using std::ostringstream;

using namespace genie::nuvld;

//______________________________________________________________________________
GuiDataSelectionDialog::GuiDataSelectionDialog(void) 
{

}
//______________________________________________________________________________
GuiDataSelectionDialog::~GuiDataSelectionDialog()
{

}
//______________________________________________________________________________
string GuiDataSelectionDialog::BundleSelectionsInString(void)
{
  ostringstream options;

  options << "KEY-LIST:" << BundleKeyListInString() << "$"
          << "CUTS:"     << BundleCutsInString()    << "$"
          << "DRAW_OPT:" << BundleDrawOptInString() << "$"
          << "DB-TYPE:vN-XSec"; // tmp

  return options.str();
}
//______________________________________________________________________________

