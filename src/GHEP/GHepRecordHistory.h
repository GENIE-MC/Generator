//____________________________________________________________________________
/*!

\class    genie::GHepRecordHistory

\brief    Holds the history of the GHEP event record as it being modified by 
          the processing steps of an event generation thread.
          The event record history can be used to step back in the generation
          sequence if a processing step is to be re-run (this the GENIE event
          generation framework equivalent of an 'Undo')

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 23, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _GHEP_RECORD_HISTORY_H_
#define _GHEP_RECORD_HISTORY_H_

#include <map>
#include <string>
#include <ostream>

using std::map;
using std::string;
using std::ostream;

namespace genie {

class GHepRecord;

class GHepRecordHistory : public map<int, GHepRecord*> {

public :

  GHepRecordHistory();
  GHepRecordHistory(const GHepRecordHistory & history);
  ~GHepRecordHistory();

  void AddSnapshot        (int step, GHepRecord * r);
  void PurgeHistory       (void);
  void PurgeRecentHistory (int start_step);
  void ReadFlags          (void);

  void Copy  (const GHepRecordHistory & history);
  void Print (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const GHepRecordHistory & history);

private:

  bool fEnabledFull;          ///< keep the full GHEP record history
  bool fEnabledBootstrapStep; ///< keep only the record that bootsrapped the generation cycle
};

}      // genie namespace

#endif // _GHEP_RECORD_HISTORY_H_
