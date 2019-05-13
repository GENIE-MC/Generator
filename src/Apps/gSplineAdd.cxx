//____________________________________________________________________________
/*!

\program gspladd

\brief   Merges XML files containing GENIE cross section splines

         Syntax :
           gspladd -f file_list -d directory_list -o output.xml
                   [--message-thresholds xml_file]

         Options :
           -f 
              A list of input xml cross-section files. If more than one then 
              separate using commas.
           -d 
              A list of input directories where to look for xml cross section
              files. If more than one then separate using commas.
           -o 
              output xml file
           --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.

         Notes :
           There must be at least 2 files for the merges to work

         Examples :

           1) shell% gspladd -f xsec_Fe56.xml,xsec_O16.xml -o xsec.xml

              will merge xsec_Fe56.xml and xsec_O16.xml into a single file
              named xsec_all.xml

           2) shell% gspladd -f xsec_Fe56.xml -d /path,/other_path  -o xsec.xml

              will merge xsec_Fe56.xml with all the xml cross-section files that
              can be found in the /path and /other_path directories and write-out
              a single file named xsec_all.xml

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         Rutherford Appleton Laboratory

\created July 05, 2007

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>
#include <vector>

#include <TSystem.h>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/XmlParserStatus.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;

vector<string> GetAllInputFiles   (void);
void           GetCommandLineArgs (int argc, char ** argv);
void           PrintSyntax        (void);

//User-specified options:
string         gOutFile;   ///< output XML file
vector<string> gInpFiles;  ///< list of input XML files
vector<string> gInpDirs;   ///< list of input dirs (to look for XML files)
vector<string> gAllFiles;  ///< list of all input files

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc,argv);

  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  
  XSecSplineList * xspl = XSecSplineList::Instance();

  vector<string>::const_iterator file_iter = gAllFiles.begin();
  for( ; file_iter != gAllFiles.end(); ++file_iter) {
    string filename = *file_iter;
    LOG("gspladd", pNOTICE) << " ---- >> Loading file : " << filename;
    XmlParserStatus_t ist = xspl->LoadFromXml(filename, true);
    assert(ist==kXmlOK);
  }

  LOG("gspladd",pDEBUG) << *xspl ;

  LOG("gspladd", pNOTICE) 
     << " ****** Saving all loaded splines into : " << gOutFile;
  xspl->SaveAsXml(gOutFile);

  return 0;
}
//____________________________________________________________________________
vector<string> GetAllInputFiles(void)
{
  vector<string> files;

  vector<string>::const_iterator file_iter;
  vector<string>::const_iterator dir_iter;

  // add all files that were input explictly
  file_iter = gInpFiles.begin();
  for( ; file_iter != gInpFiles.end(); ++file_iter) {
    string filename = *file_iter;
    files.push_back(filename);
  } // file_iter

  // loop over input directories
  dir_iter = gInpDirs.begin();
  for( ; dir_iter != gInpDirs.end(); ++dir_iter) {
    string path = *dir_iter;
    // get all XML files in this dir
    vector<string> path_files = utils::system::GetAllFilesInPath(path,"xml");
    // add these files too
    file_iter = path_files.begin();
    for( ; file_iter != path_files.end(); ++file_iter) {
       string filename = *file_iter;
       files.push_back(filename);
    }//file_iter
  }//dir_iter

  return files;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gspladd", pNOTICE) << "Parsing command line arguments";

  // Common run options. 
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  if( parser.OptionExists('f') ) {
    LOG("gspladd", pINFO) << "Reading input files";
    string inpfiles = parser.ArgAsString('f');
    if(inpfiles.find(",") != string::npos) {
       // split the comma separated list
       gInpFiles = utils::str::Split(inpfiles, ",");
    } else {
       // there is just one file
       gInpFiles.push_back(inpfiles);
    }
  }

  if( parser.OptionExists('d') ) {
    LOG("gspladd", pINFO) << "Reading input directories";
    string inpdirs = parser.ArgAsString('d');
    if(inpdirs.find(",") != string::npos) {
       // split the comma separated list
       gInpDirs = utils::str::Split(inpdirs, ",");
    } else {
       // there is just one directory
       gInpDirs.push_back(inpdirs);
    }
  }

  if( parser.OptionExists('o') ) {
    LOG("gspladd", pINFO) << "Reading output file name";
    gOutFile = parser.ArgAsString('o');
  } else {
    LOG("gspladd", pFATAL) << "You must specify an output file name";
    PrintSyntax();
    exit(1);
  }

  gAllFiles = GetAllInputFiles();
  if(gAllFiles.size() <= 1) {
    LOG("gspladd", pFATAL) << "There must be at least 2 input files";
    PrintSyntax();
    exit(1);
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gspladd", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gspladd  -f file_list -d directory_list  -o output.xml\n"
    << "            [--message-thresholds xml_file]\n";

}
//____________________________________________________________________________
