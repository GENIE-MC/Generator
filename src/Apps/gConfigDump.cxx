//____________________________________________________________________________
/*!

\program

\brief   Dump the config for a tune in human readable format


         Syntax :
           gconfigdump [-h]
                  [--event-generator-list list_name]
                  [--tune genie_tune]

         Options :
           [] Denotes an optional argument.
           -h
              Prints-out help on using gconfidump and exits.
           --event-generator-list
              List of event generators to load in event generation drivers.
              [default: "Default"].
           --tune
              Specifies a GENIE comprehensive neutrino interaction model tune.
              [default: "Default"].

        ***  See the User Manual for more details and examples. ***

\author  Robert Hatcher <rhatcher \at fnal.gov>
 Fermilab

\created September 21, 2021

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <map>

//#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
//#include <fenv.h> // for `feenableexcept`
//#endif

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TFolder.h>
#include <TObjString.h>
#include <TROOT.h>

#include "Framework/Conventions/XmlParserStatus.h"
#include "Framework/Utils/XmlParserUtils.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/EnvSnapshot.h"

#include "Framework/Messenger/Messenger.h"

#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"

#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCJobConfig.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/StringUtils.h"

using std::string;
using std::vector;
using std::map;
using std::ostringstream;

using namespace genie;
using namespace genie::controls;

void GetCommandLineArgs (int argc, char ** argv);
void Initialize         (void);
void PrintSyntax        (void);


//User-specified options:

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc,argv);
  Initialize();

  // throw on NaNs and Infs...
//#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
//#endif

/*
  Registry* rsub = new Registry;
  rsub->Set("xyz",1234.);

  RgAlg ralg("Name","Config");

  Registry* rtop = new Registry;
  rtop->Set("abc",9876.);
  //rtop->Set("rsub",rsub);
  rtop->Set("ralg",ralg);

  std::cout << *rtop << std::endl;
  return 0;
*/

  const string common_key_root = "Common";
  // key,file
  std::map<std::string,std::string> common_list;  // e.g. CKM,Param
  std::set<std::string> common_files;


  // taken from NtpMCJobConfig
  AlgConfigPool * algconf = AlgConfigPool::Instance();

  //TFolder* fConfig =
  //  gROOT->GetRootFolder()->AddFolder("gconfig","GENIE configs");

  const vector<string> & vconfkeys = algconf->ConfigKeyList();
  vector<string>::const_iterator keyiter;

  for(keyiter = vconfkeys.begin(); keyiter != vconfkeys.end(); ++keyiter) {

    string key = *keyiter;
    std::cout << "Registry key: " << key << std::endl;

    vector<string> vkey = utils::str::Split(key,"/");
    assert(vkey.size()==2);
    string alg_name  = vkey[0];
    string param_set = vkey[1];

    std::cout
      << "alg_name: " << alg_name << ", param_set: " << param_set; // << std::endl;

    //if( !(fConfig->FindObject(alg_name.c_str())) ) {
    //  LOG("gconfigdump",pDEBUG) << "Adding new folder for alg: " << alg_name;
    //  fConfig->AddFolder(alg_name.c_str(), "");
    //}
    //TFolder * alg_folder = (TFolder *) fConfig->FindObject(alg_name.c_str());

    //LOG("gconfigdump",pDEBUG) << "Adding folder for param set: " << param_set;
    //TFolder * config_folder = alg_folder->AddFolder(param_set.c_str(), "");

    //LOG("gconfigdump",pDEBUG) << "Accessing Registry & converting it to TFolder";
    Registry * config_registry = algconf->FindRegistry(key);
    //config_registry->CopyToFolder(config_folder);
    std::cout << *config_registry;

    // pick out items that start with Common
    for ( RgIMapConstIter it  = config_registry->GetItemMap().begin();
          it != config_registry->GetItemMap().end() ; ++it ) {
      std::string keynm = it->first;
      if ( keynm.find(common_key_root) == 0 ) {
        // strip "Common" off key ... leave Param or Decay
        std::string type = keynm.substr( common_key_root.size() );
        //
        std::string cpstr = config_registry->GetStringDef(keynm,"not-found",false);
        if ( cpstr != "not-found") {
          vector<std::string> cps = utils::str::Split(cpstr,",");
          for (size_t i=0; i<cps.size(); ++i) {
            common_list[cps[i]] = type;
            common_files.insert(type);
          }
        }
      }
    }

    std::cout << std::endl;
  }

  std::string xmlpaths = genie::utils::xml::GetXMLPathList();
  std::cout << "GetXMLPathList returns: "; // << xmlpaths << std::endl;
  vector<string> xmlpathlist = utils::str::Split(xmlpaths,":");
  auto xmlpathitr = xmlpathlist.begin();
  for ( ; xmlpathitr != xmlpathlist.end(); ++xmlpathitr) {
    std::cout << std::endl << "  " << *xmlpathitr;
  }
  std::cout << std::endl << std::endl;


  Registry* common;

  std::set<std::string>::iterator typeitr = common_files.begin();
  for ( ; typeitr != common_files.end(); ++typeitr ) {
    std::string fname = std::string("Common") + *typeitr + std::string(".xml");
    std::string commonpath = genie::utils::xml::GetXMLFilePath(fname);
    std::cout << "using file " << commonpath << std::endl;


    std::map<std::string,std::string>::iterator cpitr = common_list.begin();
    for ( ; cpitr != common_list.end(); ++cpitr) {
      if ( cpitr->second != *typeitr ) continue;
      std::cout << "Common" << cpitr->second << " \"" << cpitr->first << "\"";

      common = algconf->CommonList(cpitr->second,cpitr->first);
      if ( ! common ) {
        std::cout << "\n\tno Common" << cpitr->second << " \"" << cpitr->first << "\"" << std::endl;
      } else {
        std::cout << *common;
      }
    }
    std::cout << std::endl;
  }

  common = algconf->CommonList("Param","Tunable");
  if ( ! common ) {
    std::cout << "no CommonParam \"Tunable\"" << std::endl;
  } else {
    std::cout << *common;
  }

  std::cout << std::endl << std::endl;

  // taken from NtpMCJobEnv
  unsigned int ivar = 0;
  while ( kMCEnv[ivar] ) {
    ostringstream entry;
    string var   = kMCEnv[ivar];
    string value = (gSystem->Getenv(var.c_str()) ?
                    gSystem->Getenv(var.c_str()) : "UNDEFINED");

    std::cout << "$" << var << " ---> " << value << std::endl;
    ivar++;
  }


}
//____________________________________________________________________________
void Initialize()
{

  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("gconfigdump", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();

  // Initialization of random number generators, cross-section table,
  // messenger thresholds, cache file
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::CacheFile(RunOpt::Instance()->CacheFile());

}

//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gconfigdump", pINFO) << "Parsing command line arguments";

  // Common run options. Set defaults and read.
  RunOpt::Instance()->EnableBareXSecPreCalc(true);
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
      PrintSyntax();
      exit(0);
  }


  //
  // print-out the command line options
  //
  LOG("gconfigdump", pNOTICE)
     << "\n"
     << utils::print::PrintFramedMesg("gconfigdump job configuration");
  LOG("gconfigdump", pNOTICE) << "\n";

  LOG("gconfigdump", pNOTICE) << *RunOpt::Instance();

}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gconfigdump", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "\n      gconfigdump [-h]"
    << "\n              [--event-generator-list list_name]"
    << "\n              [--tune G18_02a_00_000] (or your preferred tune identifier)"
    << "\n";
}
//____________________________________________________________________________
