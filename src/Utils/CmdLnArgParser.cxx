//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jul 23, 2010 - CA
   First added in v2.7.1. Adapted from previous code to parse cmd line args.

*/
//____________________________________________________________________________

#include <cctype>
#include <cstdlib>

#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLnArgParser.h"

using namespace genie;
using namespace genie::utils;

//____________________________________________________________________________
CmdLnArgParser::CmdLnArgParser(int argc, char **argv) :
fArgc(argc),
fArgv(argv)
{

}
//____________________________________________________________________________
CmdLnArgParser::~CmdLnArgParser()
{

}
//____________________________________________________________________________
char * CmdLnArgParser::Arg(char op)
{
  bool    set       = false;
  char *  argument  = new char[1024];

  int     argc      = fArgc;
  char ** argv      = fArgv;

  while(argc>2)
  {
    LOG("CLAP", pDEBUG) << "Getting next argument in argument-list";
    LOG("CLAP", pDEBUG) << "Current argc = " << argc;

    if (argv[1][0] == '-') {      
      LOG("CLAP", pDEBUG) 
       << "Got char (argv[1][1]) following argv[1][0]='-' : " << argv[1][1];

      if (argv[1][1] == op) {
  	 LOG("CLAP", pDEBUG) << "Input option: " << op << " was matched";

         if (strlen(&argv[1][2]) ) {
            strcpy(argument,&argv[1][2]);
            set = true;
            LOG("CLAP", pINFO) 
               << "Set opt = [" << op << "] to val = [" << argument << "]";

         } else if( (argc>2) && 
                    !(argv[2][0]=='-' && isalpha(argv[2][1])) ) {

            LOG("CLAP", pDEBUG) 
              << "argc>2 and next arg not a '-' followed by an alpha char";

            argc--;
            argv++;
            strcpy(argument,&argv[1][0]);
            set = true;
            LOG("CLAP", pINFO) 
               << "Set opt = [" << op << "] to val = [" << argument << "]";
         }
      }
    }
    argc--;
    argv++;
    if(argc>2) {
      LOG("CLAP", pDEBUG) << "Next argv[1][0] = " << argv[1][0];
    }
  }

  return argument;
}
//____________________________________________________________________________
bool CmdLnArgParser::Exists(char op)
{
  bool set = false;

  int     argc = fArgc;
  char ** argv = fArgv;

  while(argc>1) {
    if(argv[1][0] == '-') {
      if (argv[1][1] == op) set = true;
    }
    argc--;
    argv++;
  }

  return set;
}
//____________________________________________________________________________
string CmdLnArgParser::ArgAsString(char op)
{
  char * argument = this->Arg(op);
  string value = string(argument);
  delete [] argument;

  return value;
}
//____________________________________________________________________________
vector<string> CmdLnArgParser::ArgAsStringTokens(char op, string delimeter)
{
  string argument = this->ArgAsString(op);

  vector<string> tokens = str::Split(argument, delimeter);
  return tokens;
}
//____________________________________________________________________________
double CmdLnArgParser::ArgAsDouble(char op)
{
  char * argument = this->Arg(op);
  double value = atof(argument);
  delete [] argument;

  return value;
}
//____________________________________________________________________________
vector<double> CmdLnArgParser::ArgAsDoubleTokens(char op, string delimeter)
{
  vector<string> strtokens = this->ArgAsStringTokens(op, delimeter);
  vector<double> tokens;
  vector<string>::const_iterator iter = strtokens.begin();
  for( ; iter != strtokens.end(); ++iter) {
    string arg = *iter;
    tokens.push_back(atof(arg.c_str()));
  }
  return tokens;
}
//____________________________________________________________________________
int CmdLnArgParser::ArgAsInt(char op)
{
  char * argument = this->Arg(op);
  int value = atoi(argument);
  delete [] argument;

  return value;
}
//____________________________________________________________________________
vector<int> CmdLnArgParser::ArgAsIntTokens(char op, string delimeter)
{
  vector<string> strtokens = this->ArgAsStringTokens(op, delimeter);
  vector<int> tokens;
  vector<string>::const_iterator iter = strtokens.begin();
  for( ; iter != strtokens.end(); ++iter) {
    string arg = *iter;
    tokens.push_back(atoi(arg.c_str()));
  }
  return tokens;
}
//____________________________________________________________________________
long CmdLnArgParser::ArgAsLong(char op)
{
  char * argument = this->Arg(op);
  long value = atol(argument);
  delete [] argument;

  return value;
}
//____________________________________________________________________________
vector<long> CmdLnArgParser::ArgAsLongTokens(char op, string delimeter)
{
  vector<string> strtokens = this->ArgAsStringTokens(op, delimeter);
  vector<long> tokens;
  vector<string>::const_iterator iter = strtokens.begin();
  for( ; iter != strtokens.end(); ++iter) {
    string arg = *iter;
    tokens.push_back(atol(arg.c_str()));
  }
  return tokens;
}
//____________________________________________________________________________
