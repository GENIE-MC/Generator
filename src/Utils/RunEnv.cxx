//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 29, 2013 - CA
   Added in preparartion for v2.8.0, where use of env. vars was phased out.

*/
//____________________________________________________________________________

//#include <sstream>
#include <iostream>

//#include <TSystem.h>

#include "Messenger/Messenger.h"
#include "Utils/RunEnv.h"

//using std::ostringstream;
using std::cout;
using std::endl;

namespace genie {

//____________________________________________________________________________
ostream & operator << (ostream & stream, const RunEnv & env)
{
  env.Print(stream);
  return stream;
}
//____________________________________________________________________________
RunEnv * RunEnv::fInstance = 0;
//____________________________________________________________________________
RunEnv::RunEnv()
{
  fInstance = 0;

  this->EnableCache(true);
}
//____________________________________________________________________________
RunEnv::~RunEnv()
{
  fInstance = 0;
}
//____________________________________________________________________________
RunEnv * RunEnv::Instance()
{
  if(fInstance == 0) {
    static RunEnv::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new RunEnv;
  }
  return fInstance;
}
//____________________________________________________________________________
void RunEnv::Print(ostream & stream) const
{
  stream << "\n Cache .... : " << ((fEnableCache) ? "Enabled" : "Disabled");
  stream << "\n";
}
//___________________________________________________________________________

} // genie namespace


