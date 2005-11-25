//____________________________________________________________________________
/*!

\class   genie::GMCJMonitor

\brief   Simple class to create & update MC job status files and env. vars.
         This is used to be able to keep track of an MC job status even when
         all output is suppressed or redirected to /dev/null.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 13, 2005

*/
//____________________________________________________________________________

#include <sstream>
#include <fstream>

#include <TSystem.h>

#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "Utils/PrintUtils.h"

using std::ostringstream;
using std::endl;
using std::ios;
using std::ofstream;

using namespace genie;

//____________________________________________________________________________
GMCJMonitor::GMCJMonitor(Long_t runnu) :
fRunNu(runnu)
{
  this->BuildNames();
}
//____________________________________________________________________________
GMCJMonitor::~GMCJMonitor()
{

}
//____________________________________________________________________________
void GMCJMonitor::Update(int iev, const EventRecord * event)
{
  this->UpdateStatusEnvVar(iev, event);
  this->UpdateStatusFile(iev, event);
}
//____________________________________________________________________________
void GMCJMonitor::UpdateStatusEnvVar(int iev, const EventRecord * event)
{
  ostringstream status;

  status << "IEv: " << iev << ";";

  if(!event) {
     status << "NULL";
  } else {
     status << "Ev:"
            << event->GetInteraction()->AsString()
            << ";";
     status << "P4:"
            << utils::print::P4AsShortString(event->GetParticle(0)->P4())
            << ";";
     status << "Vtx:"
            << utils::print::X4AsString(event->GetParticle(0)->V4())
            << ";";
  }

  LOG("GMCJMonitor", pINFO) << "$" << fStatusEnvVar << "->" << status.str();

  gSystem->Setenv(fStatusEnvVar.c_str(), status.str().c_str());
}
//____________________________________________________________________________
void GMCJMonitor::UpdateStatusFile(int iev, const EventRecord * event)
{
  ofstream out(fStatusFile.c_str(), ios::out);

  ostringstream status;

  status << endl;
  status << "Current Event Number: " << iev << endl;
  if(!event) status << "NULL" << endl;
  else       status << *event << endl;

  out << status.str();
  out.close();
}
//____________________________________________________________________________
void GMCJMonitor::BuildNames(void)
{
  ostringstream filename;
  ostringstream envvarname;

  filename   << "genie-mcjob-" << fRunNu << ".status";
  envvarname << "GMCSTATUS" << fRunNu;

  fStatusFile   = filename.str();
  fStatusEnvVar = envvarname.str();
}
//____________________________________________________________________________

