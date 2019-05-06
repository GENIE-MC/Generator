#ifndef _HINCL_test_H_
#define _HINCL_test_H_


#include <TGenPhaseSpace.h>

#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Physics/HadronTransport/INukeMode.h"
#include "Physics/HadronTransport/INukeHadroFates.h"
#include "Physics/HadronTransport/Intranuke2018.h"
#include "INCLCascade.h"

class TLorentzVector;
class TVector3;

namespace genie{
  class GHepParticle;
  class INukeHadroData;
  class PDGCodeList;
class HINCLCascade: public INCLCascade{
  friend class IntranukeTester;
public :
 HINCLCascade();
 HINCLCascade(string config);

 ~HINCLCascade();
 void ProcessEventRecord(GHepRecord * event_rec) const;
 private:
   void LoadConfig (void);
};
}
#endif
