//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - July 15, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <iostream>

#include "EVGCore/RunningThreadInfo.h"
#include "EVGCore/EventGeneratorI.h"
#include "Messenger/Messenger.h"

using std::endl;
using std::cout;

using namespace genie;

//____________________________________________________________________________
RunningThreadInfo * RunningThreadInfo::fInstance = 0;
//____________________________________________________________________________
RunningThreadInfo::RunningThreadInfo()
{
  fInstance =  0;
}
//____________________________________________________________________________
RunningThreadInfo::~RunningThreadInfo()
{
  cout << "RunningThreadInfo singleton dtor" << endl;
  fInstance = 0;
}
//____________________________________________________________________________
RunningThreadInfo * RunningThreadInfo::Instance()
{
  if(fInstance == 0) {
    static RunningThreadInfo::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new RunningThreadInfo;
  }
  return fInstance;
}
//____________________________________________________________________________
