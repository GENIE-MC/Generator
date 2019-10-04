//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 18, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TDatime.h>

#include "Framework/Ntuple/NtpMCDTime.h"

using namespace genie;

ClassImp(NtpMCDTime)

using std::endl;

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream& stream, const NtpMCDTime & dt)
  {
     dt.PrintToStream(stream);
     return stream;
  }
}
//____________________________________________________________________________
NtpMCDTime::NtpMCDTime() :
TObject()
{
  this->Init();
}
//____________________________________________________________________________
NtpMCDTime::NtpMCDTime(const NtpMCDTime & dt) :
TObject()
{
  this->Copy(dt);
}
//____________________________________________________________________________
NtpMCDTime::~NtpMCDTime()
{

}
//____________________________________________________________________________
void NtpMCDTime::PrintToStream(ostream & stream) const
{
  stream
    << "DATE (dd/mm/yyyy): "   << day  << "/" << month << "/" << year
    << ", TIME (hr:min:sec): " << hour << ":" << min   << ":" << sec << endl;
}
//____________________________________________________________________________
void NtpMCDTime::Copy(const NtpMCDTime & dt)
{
  year   = dt.year;
  month  = dt.month;
  day    = dt.day;
  hour   = dt.hour;
  min    = dt.min;
  sec    = dt.sec;
  val    = dt.val;
}
//____________________________________________________________________________
void NtpMCDTime::Init(void)
{
  year   = 0;
  month  = 0;
  day    = 0;
  hour   = 0;
  min    = 0;
  sec    = 0;
  val    = 0;
}
//____________________________________________________________________________
void NtpMCDTime::Now(void)
{
  TDatime t;

  year   = t.GetYear();
  month  = t.GetMonth();
  day    = t.GetDay();
  hour   = t.GetHour();
  min    = t.GetMinute();
  sec    = t.GetSecond();
  val    = t.Get();
}
//____________________________________________________________________________
