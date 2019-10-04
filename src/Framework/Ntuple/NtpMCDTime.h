//____________________________________________________________________________
/*!

\class   genie::NtpMCDTime

\brief   MINOS-style Ntuple Class to hold the date and time that the event
         ntuple was generated.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created October 18, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NTP_MC_DTIME_H_
#define _NTP_MC_DTIME_H_

#include <ostream>

#include <TObject.h>

using std::ostream;

namespace genie {

class NtpMCDTime;
ostream & operator << (ostream & stream, const NtpMCDTime & dt);

class NtpMCDTime : public TObject {

public :
  using TObject::Copy; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings

  NtpMCDTime();
  NtpMCDTime(const NtpMCDTime & dt);
  virtual ~NtpMCDTime();

  void Now  (void);
  void Init (void);
  void Copy (const NtpMCDTime & dt);

  void PrintToStream(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const NtpMCDTime & dt);

  // Ntuple is treated like a C-struct with public data members and
  // rule-breakinsg field data members not prefaced by "f" and mostly lowercase.
  Int_t         year;
  Int_t         month;
  Int_t         day;
  Int_t         hour;
  Int_t         min;
  Int_t         sec;
  UInt_t        val;  ///< sec since 1st Jan 1995 as in TDatime

  ClassDef(NtpMCDTime, 1)
};

}      // genie namespace

#endif // _NTP_MC_DTIME_H_
