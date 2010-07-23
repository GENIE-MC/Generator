//____________________________________________________________________________
/*!

\program gtestCmdLnArg

\brief   Program used for testing the command line argument parser

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created July 23, 2010

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>
#include <vector>

#include "Messenger/Messenger.h"
#include "Utils/CmdLnArgParser.h"

using std::string;
using std::vector;

using namespace genie;

int main(int argc, char **argv)
{
  CmdLnArgParser parser(argc,argv);

  LOG("test", pNOTICE)
     << "Option -f exists? " << parser.Exists('f');
  LOG("test", pNOTICE)
     << "Option -n exists? " << parser.Exists('n');
  LOG("test", pNOTICE)
     << "Option -s exists? " << parser.Exists('s');

  if(parser.Exists('d')) {
    LOG("test", pNOTICE)
       << "Command line argument following -d : " 
       << parser.Arg('d');
    LOG("test", pNOTICE)
       << "Command line argument following -d (as double): " 
       << parser.ArgAsDouble('d');
  }

  if(parser.Exists('l')) {
    LOG("test", pNOTICE)
       << "Command line argument following -l : " 
       << parser.Arg('l');
    LOG("test", pNOTICE)
       << "Command line argument following -l (as long int): " 
       << parser.ArgAsLong('l');
  }

  if(parser.Exists('v')) {
    LOG("test", pNOTICE)
       << "Command line argument following -v : " 
       << parser.Arg('v');
    vector<string> tokens = parser.ArgAsStringTokens('v',",");
    vector<string>::const_iterator iter = tokens.begin();
    for( ; iter != tokens.end(); ++iter) {
       LOG("test", pNOTICE) 
         << "Token: " << *iter;
    }
  }

  return 0;
}
//____________________________________________________________________________
