//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 
*/
//____________________________________________________________________________

// for exit()
#include <cstdlib>

#include <TSystem.h>

//#include "Framework/Conventions/XmlParserStatus.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/XmlParserUtils.h"

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

  // don't try to expand if no filename actually given ...
  string expandedinpfile = "";
  string fullinpfile     = "";
  if ( inpfile != "" ) {
    // expand in case of embedded env var or ~
    expandedinpfile = gSystem->ExpandPathName(inpfile.c_str());
    if (utils::system::FileExists(expandedinpfile)) {
      // use the file as given if possible
      fullinpfile   = expandedinpfile;
    } else {
      // look for file in $GXMLPATH, then $GENIE/config
      // return input name if not found any of those places (thus allowing CWD)
      fullinpfile   = genie::utils::xml::GetXMLFilePath(expandedinpfile);
    }
  }

  // file was specified & exists - load table
  if (utils::system::FileExists(fullinpfile)) {
    xspl = XSecSplineList::Instance();
    XmlParserStatus_t status = xspl->LoadFromXml(fullinpfile);
    if (status != kXmlOK) {
      LOG("AppInit", pFATAL)
         << "Problem reading file: " << expandedinpfile;
       gAbortingInErr = true;
       exit(1);
    }
  }

  // file doesn't exist
  else {
    // if one was specified, report & exit
    if (inpfile.size() > 0) {
       LOG("AppInit", pFATAL)
          << "Input cross-section file [" << inpfile << "] does not exist!\n"
          << "looked for " << expandedinpfile << " in $GXMLPATH locations ";
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

