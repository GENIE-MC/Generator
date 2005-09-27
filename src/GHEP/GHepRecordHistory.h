//____________________________________________________________________________
/*!

\class   genie::GHepRecordHistory

\brief   Holds the history of the GHEP event record within a single event
         generation sequence as a function of its processing step.
         The event record history is used to step back in the generation
         sequence if a processing step is to be re-run (the GENIE event
         generation framework equivalent of 'Undo')

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 23, 2005

*/
//____________________________________________________________________________

#ifndef _GHEP_RECORD_HISTORY_H_
#define _GHEP_RECORD_HISTORY_H_

#include <map>
#include <ostream>

using std::map;
using std::ostream;

namespace genie {

class GHepRecord;

class GHepRecordHistory : public map<int, GHepRecord*> {

public :

  GHepRecordHistory();
  GHepRecordHistory(const GHepRecordHistory & history);
  ~GHepRecordHistory();

  void AddSnapshot  (int step, GHepRecord * record);
  void PurgeHistory (void);

  void Copy  (const GHepRecordHistory & history);
  void Print (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const GHepRecordHistory & history);
};

}      // genie namespace

#endif // _GHEP_RECORD_HISTORY_H_
