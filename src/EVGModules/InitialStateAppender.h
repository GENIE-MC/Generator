//____________________________________________________________________________
/*!

\class   genie::InitialStateAppender

\brief   Appends the initial state information to the event record.

         Is a concerete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 04, 2004

*/
//____________________________________________________________________________

#ifndef _INITIAL_STATE_APPENDER_H_
#define _INITIAL_STATE_APPENDER_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class InitialStateAppender : public EventRecordVisitorI {

public :

  InitialStateAppender();
  InitialStateAppender(const char * param_set);
  ~InitialStateAppender();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

private :

  void AddNeutrino       (GHepRecord * event_rec) const;
  void AddNucleus        (GHepRecord * event_rec) const;
  void AddStruckParticle (GHepRecord * event_rec) const;

};

}      // genie namespace

#endif // _INITIAL_STATE_APPENDER_H_
