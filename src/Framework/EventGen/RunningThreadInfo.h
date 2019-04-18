//____________________________________________________________________________
/*!

\class    genie::RunningThreadInfo

\brief    Keep info on the event generation thread currently on charge.
          This is used so that event generation modules invoked by the thread
	  can see the "bigger picture" and access the cross section model for
	  the thread, look-up info for modules that run before or are scheduled
          to run after etc.
	  
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 06, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RUNNING_THREAD_INFO_H_
#define _RUNNING_THREAD_INFO_H_

namespace genie {

class EventGeneratorI;

class RunningThreadInfo
{
public:
  static RunningThreadInfo * Instance(void);

  const EventGeneratorI * RunningThread(void) 
  {   
    return fRunningThread; 
  }
  void UpdateRunningThread(const EventGeneratorI * evg) 
  { 
     fRunningThread = evg; 
  }

private:
  RunningThreadInfo();
  RunningThreadInfo(const RunningThreadInfo & info);
  virtual ~RunningThreadInfo();

  //! self
  static RunningThreadInfo * fInstance;

  //! current thread
  const EventGeneratorI * fRunningThread;

  //! clean
  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (RunningThreadInfo::fInstance !=0) {
            delete RunningThreadInfo::fInstance;
            RunningThreadInfo::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace

#endif // _RUNNING_THREAD_INFO_H_
