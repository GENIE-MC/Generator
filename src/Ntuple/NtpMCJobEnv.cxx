//____________________________________________________________________________
/*!

\class   genie::NtpMCJobEnv

\brief   Stores a snapshot of your environment in ROOT TFolder along with the
         output event tree

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include <string>
#include <sstream>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>
#include <TFolder.h>
#include <TObjString.h>

#include "Conventions/Controls.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCJobEnv.h"

using std::string;
using std::ostringstream;
using std::vector;

using namespace genie;
using namespace genie::controls;

//____________________________________________________________________________
NtpMCJobEnv::NtpMCJobEnv()
{
  fEnv = 0;
}
//____________________________________________________________________________
NtpMCJobEnv::~NtpMCJobEnv()
{

}
//____________________________________________________________________________
TFolder * NtpMCJobEnv::TakeSnapshot(void)
{
  if (fEnv) delete fEnv;
  fEnv = 0;

  LOG("NtpMCEnv", pNOTICE)
                << "Taking environment snapshot and saving it in a TFolder";

  fEnv = gROOT->GetRootFolder()->AddFolder("genv","GENIE user environment");
  gROOT->GetListOfBrowsables()->Add(fEnv,"genv");

  fEnv->SetOwner(true);

  vector<string> envvars;

  unsigned int ivar=0;
  while(kMCEnv[ivar]) {
    envvars.push_back(string(kMCEnv[ivar]));
    ivar++;
  }

  LOG("NtpMCEnv", pINFO) << "** MC Job Environment:";
  vector<string>::iterator variter;
  for(variter = envvars.begin(); variter != envvars.end(); ++variter) {

     ostringstream entry;
     string var   = *variter;
     string value = (gSystem->Getenv(var.c_str()) ?
                               gSystem->Getenv(var.c_str()) : "UNDEFINED");

     LOG("NtpMCEnv", pINFO) << "$" << var << " ---> " << value;

     entry << "envv:" << var << ";value:" << value;
     fEnv->Add(new TObjString(entry.str().c_str()));
  }
  return fEnv;
}
//____________________________________________________________________________
