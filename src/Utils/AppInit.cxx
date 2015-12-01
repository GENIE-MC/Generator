//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab -

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 31, 2013 - CA
   Added in preparation for v2.8.0

*/
//____________________________________________________________________________

// for exit()
#include <cstdlib>

#include <TSystem.h>

//#include "Conventions/XmlParserStatus.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/Cache.h"
#include "Utils/XSecSplineList.h"
#include "Utils/SystemUtils.h"
#include "Utils/AppInit.h"
#include "Utils/StringUtils.h"

using namespace genie;

//___________________________________________________________________________
void genie::utils::app_init::RandGen(long int seed)
{
  // Set random number seed, if a value was set at the command-line.
  if(seed > 0) {
    RandomGen::Instance()->SetSeed(seed);
  }
}
//___________________________________________________________________________
void genie::utils::app_init::XSecTable (string inpfile, bool require_table)
{
  // Load cross-section splines using file specified at the command-line.

  XSecSplineList * xspl = XSecSplineList::Instance();
  xspl->AutoLoad(); // display warning for usage $GSPLOAD no longer supported

  // expand in case of embedded env var or ~
  string expandedinpfile = gSystem->ExpandPathName(inpfile.c_str());

  // file was specified & exists - load table
  if(utils::system::FileExists(expandedinpfile)) {
    XSecSplineList * xspl = XSecSplineList::Instance();
    XmlParserStatus_t status = xspl->LoadFromXml(expandedinpfile);
    if(status != kXmlOK) {
      LOG("AppInit", pFATAL)
         << "Problem reading file: " << expandedinpfile;
       gAbortingInErr = true;
       exit(1);
    }
  } 

  // file doesn't exist
  else {
    // if one was specified, report & exit
    if(inpfile.size() > 0) {
       LOG("AppInit", pFATAL)
          << "Input cross-section file [" << inpfile << "] does not exist!";
       gAbortingInErr = true;
       exit(1);
    } 
    // if one was not specified, warn and decide whether to exit based on the
    // input `require_table' flag
    else {
      if(!require_table) {
         LOG("AppInit", pWARN) << "No cross-section file was specified in the application inputs";
         LOG("AppInit", pWARN) << "If none is loaded, event generation might be inefficient";
      } else {
         LOG("AppInit", pFATAL) << "No cross-section file was specified in the application inputs";
         LOG("AppInit", pFATAL) << "This is mandatory as, otherwise, event generation will be prohibitively inefficient";
         gAbortingInErr = true;
         exit(1);
      }
    }
  }

}
//___________________________________________________________________________
void genie::utils::app_init::MesgThresholds(string filelist)
{
  std::vector<std::string> files = genie::utils::str::Split(filelist,":;,");
  for (size_t i=0; i < files.size(); ++i ) {
    std::string inp_file = files[i];
    if(inp_file.size() > 0) {
      Messenger * m = Messenger::Instance();
      bool ok = m->SetPrioritiesFromXmlFile(inp_file);
      if(!ok) {
        LOG("AppInit", pWARN) 
          << "Could not load customized mesg thresholds from: " 
          << inp_file;
      }
    }
  }

}
//___________________________________________________________________________
void genie::utils::app_init::CacheFile(string inp_file)
{
  if(inp_file.size() > 0) {
    Cache::Instance()->OpenCacheFile(inp_file);
  }
}
//___________________________________________________________________________

