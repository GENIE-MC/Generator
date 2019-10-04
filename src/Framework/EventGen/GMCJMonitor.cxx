//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jun 23, 2008 - CA
   fCpuTime wasn't initialized.
 @ Jan 30, 2013 - CA
   Added SetRefreshRate(int rate)

*/
//____________________________________________________________________________

#include <sstream>
#include <fstream>
#include <cstdlib>

#include <TSystem.h>
#include <TMath.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/GMCJMonitor.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/PrintUtils.h"

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
void GMCJMonitor::SetRefreshRate(int rate) 
{ 
  fRefreshRate = TMath::Max(1,rate); 
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
  fWatch.Reset(); 
  fWatch.Start();
  fCpuTime = 0;

  // get rehreah rate of set default / protect from invalid refresh rates
  if( gSystem->Getenv("GMCJMONREFRESH") ) {
   fRefreshRate = atoi( gSystem->Getenv("GMCJMONREFRESH") );
  } else fRefreshRate = 100;

  fRefreshRate = TMath::Max(1,fRefreshRate);
}
//____________________________________________________________________________

void GMCJMonitor::CustomizeFilename(string filename)
{
  fStatusFile = filename;
}
//____________________________________________________________________________
