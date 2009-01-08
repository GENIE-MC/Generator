//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiHelpHandler

\brief    Responds to GUI events associated with the help menu

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _GUI_HELP_HANDLER_H_
#define _GUI_HELP_HANDLER_H_

#include <string>

class TGWindow;

using std::string;

namespace genie {
namespace nuvld {

class GuiHelpHandler {

public:

   GuiHelpHandler();
   GuiHelpHandler(const TGWindow * main);
   ~GuiHelpHandler();

   void NuVldAbout     (void); 
   void NuVldOnline    (void); 
   void DurhamOnline   (void); 
   void HowtoFillDBase (void); 
   void HowtoConnDBase (void);
   void Howto          (string filename); 

private:

   const TGWindow * _main;
};

} // nuvld namespace
} // genie namespace

#endif

