#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef _INCLCascade_H_
#define _INCLCascade_H_

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
  class INCL;
  class IDeExcitation;
}

namespace genie {

  class GHepParticle;
  class INukeHadroData;
  class PDGCodeList;
  class HINCLCascade;

  class INCLCascade: public EventRecordVisitorI {

  public :
    INCLCascade();
    INCLCascade(string name);
    INCLCascade(string name, string config);
    ~INCLCascade();

    int doCascade(int nflags, const char * flags[], GHepRecord * event_rec) const;

    int pdgcpiontoA(int pdgc) const;
    int pdgcpiontoZ(int pdgc) const;

    // implement the EventRecordVisitorI interface

    void Configure (const Registry & config);
    void Configure (string param_set);

    virtual void ProcessEventRecord(GHepRecord * event_rec) const;

  protected:
    bool CanRescatter(const GHepParticle * p) const;
    bool IsInNucleus(const GHepParticle * p) const;
    void TransportHadrons(GHepRecord * evrec) const;
    bool NeedsRescattering(const GHepParticle * p) const;
    virtual void LoadConfig (void) = 0;

    mutable int            fRemnA;         ///< remnant nucleus A
    mutable int            fRemnZ;         ///< remnant nucleus Z
    mutable double         fTrackingRadius;
    mutable TLorentzVector fRemnP4;        ///< P4 of remnant system
    mutable GEvGenMode_t   fGMode;
    double       fR0;           ///< effective nuclear size param
    double       fNR;

    mutable G4INCL::INCL *theINCLModel;
    mutable G4INCL::IDeExcitation *theDeExcitation;

  };

}

#endif
#endif // __GENIE_INCL_ENABLED__
