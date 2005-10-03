//____________________________________________________________________________
/*!

\namespace  genie::utils::clap

\brief      Simple utility methods for Command Line Argument Parsing in GENIE
            standalone applications

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    October 3, 2005

*/
//____________________________________________________________________________

#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"
#include "Messenger/Messenger.h"

//____________________________________________________________________________
char * genie::utils::clap::CmdLineArg(int argc, char ** argv, char op)
{
  bool   set      = false;
  char * argument = new char[1024];

  char ** argv_keep = argv;

  while(argc>2)
  {
    LOG("CLAP", pDEBUG) << "Getting next argument in argument-list";

    if (argv[1][0] == '-') {
      LOG("CLAP", pDEBUG) << "argv[1][1] = " << argv[1][1];
      if (argv[1][1] == op) {
         if (strlen(&argv[1][2]) ) {
            strcpy(argument,&argv[1][2]);
            set = true;
            LOG("CLAP", pINFO) << "Set opt = [" << op << "] "
                                  "to value = [" << argument << "]";
         } else if( (argc>2) && (argv[2][0] != '-') ) {
            argc--;
            argv++;
            strcpy(argument,&argv[1][0]);
            set = true;
            LOG("CLAP", pINFO) << "Set opt = [" << op << "] "
                                  "to value = [" << argument << "]";
         }
      }
    }
    argc--;
    argv++;
    LOG("CLAP", pDEBUG) << "Current argc    = " << argc;
    LOG("CLAP", pDEBUG) << "Next argv[1][0] = " << argv[1][0];
  }

  LOG("CLAP", pINFO) << "Restoring argv pointer";
  argv = argv_keep;

  if(!set) {
     delete [] argument;
     LOG("CLAP", pINFO)
            << "Option [" << op << "] not set. "
                            << "Throwing a CmdLineArgParserException";
     CmdLineArgParserException exception;
     exception.SetReason("Command line argument not found");
     exception.SwitchOnArgNotFound();
     throw exception;
     return 0;
  }
  return argument;
}
//____________________________________________________________________________
string genie::utils::clap::CmdLineArgAsString(int argc, char ** argv, char op)
{
  char * argument = CmdLineArg(argc, argv, op);
  string arg_as_string = string(argument);
  delete [] argument;

  return arg_as_string;
}
//____________________________________________________________________________
double genie::utils::clap::CmdLineArgAsDouble(int argc, char ** argv, char op)
{
  char * argument = CmdLineArg(argc, argv, op);
  double arg_as_double = atof(argument);
  delete [] argument;

  return arg_as_double;
}
//____________________________________________________________________________
int genie::utils::clap::CmdLineArgAsInt(int argc, char ** argv, char op)
{
  char * argument = CmdLineArg(argc, argv, op);
  int arg_as_int = atoi(argument);
  delete [] argument;

  return arg_as_int;
}
//____________________________________________________________________________
bool genie::utils::clap::CmdLineArgAsBool(int argc, char ** argv, char op)
{
  bool set = false;
  char ** argv_keep = argv;

  while(argc>1) {
    if(argv[1][0] == '-') {
      if (argv[1][1] == op) set = true;
    }
    argc--;
    argv++;
  }

  LOG("CLAP", pINFO) << "Restoring argv pointer";
  argv = argv_keep;

  return set;
}
//____________________________________________________________________________
