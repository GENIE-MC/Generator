//____________________________________________________________________________
/*!

\program gspladd

\brief   Merges XML files containing GENIE cross section splines

         Syntax :
           gspladd -i input1.xml,input2.xml,... -o output.xml

         Options :
           -i a comma separated list of input xml files (can have an arbitrary
              number of files)
           -o output xml file
 
         Example:
           gspladd -i xsec_numu_Fe56.xml,xsec_numu_O16.xml -o xsec_all.xml

           will merge xsec_numu_Fe56.xml and xsec_numu_O16.xml into a single
           file named xsec_all.xml

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         Rutherford Appleton Laboratory

\created July 05, 2007

\cpright Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>
#include <vector>

#include <TSystem.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/XmlParserStatus.h"
#include "Messenger/Messenger.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

//User-specified options:
vector<string> gInpFiles;
string         gOutFile;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc,argv);
  
  XSecSplineList * xspl = XSecSplineList::Instance();

  vector<string>::const_iterator file_iter = gInpFiles.begin();
  for( ; file_iter != gInpFiles.end(); ++file_iter) {
    string filename = *file_iter;
    LOG("gspladd", pNOTICE) << " ****** Loading file : " << filename;
    XmlParserStatus_t  ist = xspl->LoadFromXml(filename, true);
    assert(ist==kXmlOK);
  }

  LOG("gspladd", pNOTICE) 
          << " ****** Saving all loaded splines into : " << gOutFile;
  xspl->SaveAsXml(gOutFile);

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gspladd", pNOTICE) << "Parsing command line arguments";

  try {
    LOG("gspladd", pINFO) << "Reading input files";
    string inpfiles = genie::utils::clap::CmdLineArgAsString(argc,argv,'i');
    if(inpfiles.find(",") != string::npos) {
       // split the comma separated list
       gInpFiles = utils::str::Split(inpfiles, ",");
       assert(gInpFiles.size()>1);
    } else {
      LOG("gspladd", pFATAL) 
         << "You should specify at least 2 input files for the file merger to work...";
      PrintSyntax();
      exit(1);
    }
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gspladd", pFATAL) << "You must specify a list of input files";
      PrintSyntax();
      exit(1);
    }
  }


  try {
    LOG("gspladd", pINFO) << "Reading output file name";
    gOutFile = genie::utils::clap::CmdLineArgAsString(argc,argv,'o');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gspladd", pFATAL) << "You must specify an output file name";
      PrintSyntax();
      exit(1);
    }
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gspladd", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gspladd -i input1.xml,impu2.xml,... -o output.xml\n";
}
//____________________________________________________________________________
