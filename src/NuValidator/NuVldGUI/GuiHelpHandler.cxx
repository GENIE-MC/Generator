//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiHelpHandler

\brief    Responds to GUI events associated with the help menu

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include <vector>

#include <TSystem.h>
#include <TGWindow.h>

#include "NuVldGUI/GuiHelpHandler.h"
#include "NuVldGUI/HelpBox.h"
#include "NuVldGUI/SysLogSingleton.h"
#include "NuVldGUI/MultiLineMsgBox.h"
#include "Utils/StringUtils.h"

using std::vector;

using namespace genie::string_utils;
using namespace genie::nuvld;

//______________________________________________________________________________
GuiHelpHandler::GuiHelpHandler()
{
  _main = 0;
}
//______________________________________________________________________________
GuiHelpHandler::GuiHelpHandler(const TGWindow * main) : _main(main)
{

}
//______________________________________________________________________________
GuiHelpHandler::~GuiHelpHandler()
{

}
//______________________________________________________________________________
void GuiHelpHandler::NuVldAbout(void)
{
  SysLogSingleton * syslog = SysLogSingleton::Instance();

  syslog -> Log()       -> AddLine( "Opening Help->About" );
  syslog -> StatusBar() -> SetText( "Opening Help->About", 0 );

  vector<string> about;

  about.push_back("                                                         ");
  about.push_back("                 NuValidator GUI Prototype               ");
  about.push_back("                                                         ");
  about.push_back("   Costas Andreopoulos (RAL) and Hugh Gallagher (Tufts)  ");
  about.push_back("                                                         ");

  new MultiLineMsgBox(gClient->GetRoot(), _main,
                                             380, 250, kVerticalFrame, &about);

}
//______________________________________________________________________________
void GuiHelpHandler::NuVldOnline(void)
{
  int          status = -1;
  const char * url =
             "http://hepunx.rl.ac.uk/~candreop/generators/validator/index.html";

  status = gSystem->Exec( Concat("mozilla ", url) );

  if(status != 0) {
       status = gSystem->Exec( Concat("netscape ", url) );
  }
}
//______________________________________________________________________________
void GuiHelpHandler::DurhamOnline(void)
{
  gSystem->Exec(
       "mozilla http://h2.phyip3.dur.ac.uk/hepdata/online/neutrino/index.html");
}
//______________________________________________________________________________
void GuiHelpHandler::HowtoFillDBase(void)
{
  string filename = "fdb.hlp";

  this->Howto(filename);
}
//______________________________________________________________________________
void GuiHelpHandler::HowtoConnDBase(void)
{
  string filename = "ldb.hlp";

  this->Howto(filename);
}
//______________________________________________________________________________
void GuiHelpHandler::Howto(string filename)
{
  string full_filename =  string(gSystem->Getenv("GENIE")) +
                                    "/src/NuValidator/NuVldGUI/hlp/" + filename;

  new HelpBox( gClient->GetRoot(), _main, 
                              450, 200, kVerticalFrame, full_filename.c_str() );
}
//______________________________________________________________________________

