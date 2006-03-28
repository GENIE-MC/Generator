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
#include <cstdlib>

#include <TSystem.h>
#include <TMath.h>

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
  this->Init();
}
//____________________________________________________________________________
GMCJMonitor::~GMCJMonitor()
{

}
//____________________________________________________________________________
void GMCJMonitor::Update(int iev, const EventRecord * event)
{
  if(iev%fRefreshRate) return; // continue only every fRefreshRate events 

  fWatch.Stop();
  fCpuTime += (fWatch.CpuTime());

  ofstream out(fStatusFile.c_str(), ios::out);

  ostringstream status;

  status << endl;
  status << "Current Event Number: " << iev << endl;

  status << "Approximate total processing time: " 
                             << fCpuTime << " s" << endl;
  status << "Approximate processing time/event: " 
                     << fCpuTime/(iev+1) << " s" << endl;

  if(!event) status << "NULL" << endl;
  else       status << *event << endl;

  out << status.str();
  out.close();

  fWatch.Start();
}
//____________________________________________________________________________
void GMCJMonitor::Init(void)
{
  // build the filename of the GENIE status file
  ostringstream filename;
  filename   << "genie-mcjob-" << fRunNu << ".status";
  fStatusFile = filename.str();

  // create a stopwatch
  TStopwatch fWatch;
  fWatch.Reset(); 
  fWatch.Start();

  // get rehreah rate of set default / protect from invalid refresh rates
  if( gSystem->Getenv("GMCJMONREFRESH") ) {
   fRefreshRate = atoi( gSystem->Getenv("GMCJMONREFRESH") );
  } else fRefreshRate = 100;

  fRefreshRate = TMath::Max(1,fRefreshRate);
}
//____________________________________________________________________________

