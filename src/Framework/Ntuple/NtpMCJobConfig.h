//____________________________________________________________________________
/*!

\class   genie::NtpMCJobConfig

\brief   Stores the GENIE configuration in ROOT TFolders along with the
         output event tree

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created October 1, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org    
*/
//____________________________________________________________________________

#ifndef _NTP_MC_JOB_CONFIG_H_
#define _NTP_MC_JOB_CONFIG_H_

class TFolder;

namespace genie {

class NtpMCJobConfig {

public :

  NtpMCJobConfig();
  virtual ~NtpMCJobConfig();

  TFolder * Load      (void);
  TFolder * GetFolder (void) { return fConfig; }

private:

  TFolder * fConfig;
};

}      // genie namespace

#endif // _NTP_MC_JOB_CONFIG_H_
