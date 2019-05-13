//____________________________________________________________________________
/*!

\class   genie::GMCJMonitor

\brief   Simple class to create & update MC job status files and env. vars.
         This is used to be able to keep track of an MC job status even when
         all output is suppressed or redirected to /dev/null.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created July 13, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_MC_JOB_MONITOR_H_
#define _G_MC_JOB_MONITOR_H_

#include <TStopwatch.h>

namespace genie {

class EventRecord;

class GMCJMonitor {

public :
  GMCJMonitor(Long_t runnu);
 ~GMCJMonitor();

  void SetRefreshRate (int rate);
  void Update (int iev, const EventRecord * event);
  void CustomizeFilename(string filename);

private:

  void Init (void);

  Long_t     fRunNu;       ///< run number
  string     fStatusFile;  ///< name of output status file
  TStopwatch fWatch;       
  double     fCpuTime;     ///< total cpu time so far
  int        fRefreshRate; ///< update output every so many events
};

}      // genie namespace

#endif // _G_MC_JOB_MONITOR_H_
