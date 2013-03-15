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

#include <vector>

#include <TSystem.h>
#include <TGWindow.h>

#include "ValidationTools/NuVld/GuiHelpHandler.h"
#include "ValidationTools/NuVld/GuiHelpBox.h"
#include "ValidationTools/NuVld/GuiSysLogSingleton.h"
#include "ValidationTools/NuVld/GuiMultiLineMsgBox.h"
#include "Utils/StringUtils.h"

using std::vector;

using namespace genie::utils::str;
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
  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

  syslog -> Log()       -> AddLine( "Opening Help->About" );
  syslog -> StatusBar() -> SetText( "Opening Help->About", 0 );

  vector<string> about;

  about.push_back("                                                         ");
  about.push_back("                 NuValidator GUI Prototype               ");
  about.push_back("                                                         ");
  about.push_back("   Costas Andreopoulos (RAL) and Hugh Gallagher (Tufts)  ");
  about.push_back("                                                         ");

  new GuiMultiLineMsgBox(
         gClient->GetRoot(), _main, 380, 250, kVerticalFrame, &about);
}
//______________________________________________________________________________
void GuiHelpHandler::NuVldOnline(void)
{
  int status = -1;
  const char * url = "http://www.genie-mc.org";

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
  string full_filename =  
      string(gSystem->Getenv("GENIE")) +
      "/src/ValidationTools/NuVld/hlp/" + filename;

  new GuiHelpBox( 
      gClient->GetRoot(), _main, 
      450, 200, kVerticalFrame, full_filename.c_str() );
}
//______________________________________________________________________________

