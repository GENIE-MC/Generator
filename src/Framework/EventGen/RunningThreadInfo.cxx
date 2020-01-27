//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory 
*/
//____________________________________________________________________________

#include "Framework/EventGen//RunningThreadInfo.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/Messenger/Messenger.h"

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
