//____________________________________________________________________________
/*!

\class   genie::NtpMCJobEnv

\brief   Stores a snapshot of your environment in ROOT TFolder along with the
         output event tree

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created September 10, 2006

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NTP_MC_JOB_ENV_H_
#define _NTP_MC_JOB_ENV_H_

class TFolder;

namespace genie {

class NtpMCJobEnv {

public :

  NtpMCJobEnv();
  virtual ~NtpMCJobEnv();

  TFolder * TakeSnapshot (void);
  TFolder * GetFolder    (void) { return fEnv; }

private:

  TFolder * fEnv;
};

}      // genie namespace

#endif // _NTP_MC_JOB_ENV_H_
