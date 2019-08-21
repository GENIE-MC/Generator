#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef _HINCLCascadeIntranuke_H_
#define _HINCLCascadeIntranuke_H_

#include <string>
#include <TGenPhaseSpace.h>
#include <list>
#include <sstream>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>

#include "Physics/NuclearState/NuclearModelI.h"

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Conventions/GMode.h"
#include "Physics/HadronTransport/INukeMode.h"
#include "Physics/HadronTransport/INukeHadroFates.h"

class TLorentzVector;
class TVector3;
//string kDefOptEvFilePrefix="gntp.inuke";

//using namespace G4INCL;
namespace G4INCL {
  class Config;
  class INCL;
  class IDeExcitation;
}

namespace genie {

  class GHepParticle;
  class INukeHadroData;
  class PDGCodeList;

  class HINCLCascadeIntranuke: public EventRecordVisitorI {

  public :
    HINCLCascadeIntranuke();
    HINCLCascadeIntranuke(string config);
    ~HINCLCascadeIntranuke();

    int pdgcpiontoA(int pdgc) const;
    int pdgcpiontoZ(int pdgc) const;

    // implement the EventRecordVisitorI interface
    // also the LoadConfig interface

    void Configure (const Registry & config);
    void Configure (string param_set);

    virtual void ProcessEventRecord(GHepRecord * event_rec) const;

  protected:
    virtual void LoadConfig (void);

    bool CanRescatter(const GHepParticle * p) const;
    bool IsInNucleus(const GHepParticle * p) const;
    void TransportHadrons(GHepRecord * evrec) const;
    int  doCascade(GHepRecord * event_rec) const;
    bool NeedsRescattering(const GHepParticle * p) const;

    bool AddDataPathFlags(size_t& nflags, char** flags);

    mutable int            fRemnA;         ///< remnant nucleus A
    mutable int            fRemnZ;         ///< remnant nucleus Z
    mutable TLorentzVector fRemnP4;        ///< P4 of remnant system
    mutable GEvGenMode_t   fGMode;

    mutable G4INCL::Config        *theINCLConfig;
    mutable G4INCL::INCL          *theINCLModel;
    mutable G4INCL::IDeExcitation *theDeExcitation;

  };

}

#endif
#endif // __GENIE_INCL_ENABLED__
