//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, Rutherford Lab.
         October 07, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>
#include <string>
#include <vector>

#include "Viewer/GHepPrinter.h"
#include "Messenger/Messenger.h"

using std::ostringstream;
using std::string;
using std::vector;

using namespace genie;

//______________________________________________________________________________
GHepPrinter::GHepPrinter() 
{
  fGHep = 0;
}
//______________________________________________________________________________
GHepPrinter::~GHepPrinter()
{

}
//______________________________________________________________________________
void GHepPrinter::Print(EventRecord * ev_rec) 
{
  if(fGHep) {
    
    ostringstream ghep;
  
    ghep << *ev_rec;

    string ghepstr = ghep.str(); // GHEP record as a single string

    //-- Split STDHEP string to get 1 line per particle - use '\n' as delimiter

    vector<string> lines;
    string delim = string("\n");
  
    while(ghepstr.find_first_of(delim) < ghepstr.length()) {

       lines.push_back(ghepstr.substr(0, ghepstr.find_first_of(delim)) );
                          
       ghepstr = ghepstr.substr(
                       ghepstr.find_first_of(delim)+1, ghepstr.length());
     }
     lines.push_back(ghepstr);

     //-- Print GHep entries to TGTextView
   
     vector<string>::iterator line_iter;
     
     for(line_iter = lines.begin(); line_iter != lines.end(); ++line_iter) {
       fGHep->AddLine( line_iter->c_str() );
     }
  }     
}
//______________________________________________________________________________
