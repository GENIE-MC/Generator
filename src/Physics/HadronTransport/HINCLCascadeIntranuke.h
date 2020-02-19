#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef _HINCLCascadeIntranuke_H_
#define _HINCLCascadeIntranuke_H_

#include <string>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Conventions/GMode.h"

#include <TLorentzVector.h>

namespace G4INCL {
  class Config;
  class INCL;
  class IDeExcitation;
}

namespace genie {

  class GHepParticle;

  class HINCLCascadeIntranuke: public EventRecordVisitorI {

  public :
    HINCLCascadeIntranuke();
    HINCLCascadeIntranuke(std::string config);
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
    bool LookForAndAddValidPath(std::vector<std::string>& datapaths,
                                size_t defaultIndx,
                                const char* optString,
                                size_t& nflags, char** flags);


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
