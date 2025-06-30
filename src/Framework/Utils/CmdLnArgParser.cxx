//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool 
*/
//____________________________________________________________________________

#include <cctype>
#include <cstdlib>

#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"

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
  const int buf_size = 2048*128;

  bool set = false;
  char *  argument  = new char[buf_size];
  strcpy(argument, "");

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

  LOG("CLAP", pDEBUG) << "CmdLnArgParser::Arg op='" << op << "' set=" << set;
  return argument;
}
//____________________________________________________________________________
bool CmdLnArgParser::OptionExists(char op)
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
char * CmdLnArgParser::Arg(string op)
{
  const int buf_size = 2048*128;
  char *  argument   = new char[buf_size];
  strcpy(argument, "");

  int     argc      = fArgc;
  char ** argv      = fArgv;

  while(argc>2)
  {
    LOG("CLAP", pDEBUG) << "Getting next argument in argument-list";
    LOG("CLAP", pDEBUG) << "Current argc = " << argc;

    if (argv[1][0] == '-' && argv[1][1] == '-') {
      //char * op_cur = strndup(argv[1]+2,strlen(argv[1]));
      char op_cur[buf_size];
      strcpy(op_cur,&argv[1][2]);
      LOG("CLAP", pDEBUG)
        << "Got string following '--' : " << op_cur;
      if (strcmp(op.c_str(),op_cur)==0) {
  	 LOG("CLAP", pDEBUG) << "Input option: " << op << " was matched";

         if (strlen(&argv[2][0]) ) {
            strcpy(argument,&argv[2][0]);
            LOG("CLAP", pINFO)
               << "Set opt = [" << op << "] to val = [" << argument << "]";

         } else if( (argc>2) &&
                    !(argv[2][0]=='-' &&argv[2][1]=='-' && isalpha(argv[2][2])) ) {
            LOG("CLAP", pDEBUG)
              << "argc>2 and next arg not a '--' followed by an alpha char";

            argc--;
            argv++;
            strcpy(argument,&argv[1][0]);
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
bool CmdLnArgParser::OptionExists(string op)
{
  const int buf_size = 2048*128;
  bool set = false;

  int     argc = fArgc;
  char ** argv = fArgv;

  while(argc>1) {
    if(argv[1][0] == '-' && argv[1][1] == '-') {
        //char * op_cur =  strndup(argv[1]+2,strlen(argv[1]));
        char op_cur[buf_size];
        strcpy(op_cur,&argv[1][2]);
        if (strcmp(op.c_str(),op_cur)==0) set = true;
    }
    argc--;
    argv++;
  }

  return set;
}
//____________________________________________________________________________
string CmdLnArgParser::ArgAsString(string op)
{
  char * argument = this->Arg(op);
  string value = string(argument);
  delete [] argument;

  return value;
}
//____________________________________________________________________________
double CmdLnArgParser::ArgAsDouble(string op)
{
  char * argument = this->Arg(op);
  double value = atof(argument);
  delete [] argument;

  return value;
}
//____________________________________________________________________________
int CmdLnArgParser::ArgAsInt(string op)
{
  char * argument = this->Arg(op);
  int value = atoi(argument);
  delete [] argument;

  return value;
}
//____________________________________________________________________________
long CmdLnArgParser::ArgAsLong(string op)
{
  char * argument = this->Arg(op);
  long value = atol(argument);
  delete [] argument;

  return value;
}
//____________________________________________________________________________
