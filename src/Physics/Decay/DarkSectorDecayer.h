//____________________________________________________________________________
/*!
\class    genie::DarkSectorDecayer
\brief    Dark Sector decayer module.

          A simple decay simulation...
          ....
          Is a concerete implementation of the EventRecordVisitorI interface.

\author   Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
          University of Sussex

          Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  July XX, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _DARK_SECTOR_DECAYER_H_
#define _DARK_SECTOR_DECAYER_H_

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/GHEP/GHepStatus.h"

using DecayChannel = std::pair<std::vector<int>, double>;
// first the vector of pdgs, second the decay amplitude


namespace genie {

  class GHepParticle;
  class DarkSectorDecayer : public EventRecordVisitorI {

  public:
    DarkSectorDecayer();
    DarkSectorDecayer(string config);
    virtual ~DarkSectorDecayer();
    virtual void Configure(const Registry & config);
    virtual void Configure(string config);

    // Implement the EventRecordVisitorI interface
    void ProcessEventRecord(GHepRecord * event) const;

  protected:
    virtual void LoadConfig(void);

  private:
    bool ToBeDecayed(const GHepParticle & p) const;
    // comments later
    std::vector<GHepParticle> Decay(const GHepParticle & mother,
                                    const std::vector<int> & pdg_daughters) const;
    // this function will take care of the momentum conservation
    // the output particles cannot be inserted in the event record as they are
    // they need to be translated in space and time, according to the decay amplitude
    int SelectDecayChannel(const GHepParticle & mother,
                           const GHepRecord * event,
                           const std::vector<DecayChannel> & dcs,
                           double total_amplitude) const;
    std::vector<DecayChannel> DarkMediatorDecayChannels(const GHepParticle & mother,
                                                        const GHepRecord * event) const;
    std::vector<DecayChannel> DarkNeutrinoDecayChannels(const GHepParticle & mother,
                                                        const GHepRecord * event) const;
    void SetSpaceTime(std::vector<GHepParticle> & pp, const GHepParticle & mother,
                      double amplitude) const;

    mutable TGenPhaseSpace fPhaseSpaceGenerator;
    mutable double         fWeight;

    double fEps2;
    double fTheta2;
    std::array<double, 4> fMixing2s;
    double fgD2;

    double fDNuMass, fDNuMass2;
    double fDMediatorMass, fDMediatorMass2;

  };

}         // genie namespace
#endif    // _DARK_SECTOR_DECAYER_H_
