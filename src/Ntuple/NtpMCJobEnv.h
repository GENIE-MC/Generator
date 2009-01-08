//____________________________________________________________________________
/*!

\class   genie::NtpMCJobEnv

\brief   Stores a snapshot of your environment in ROOT TFolder along with the
         output event tree

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created September 10, 2006

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
