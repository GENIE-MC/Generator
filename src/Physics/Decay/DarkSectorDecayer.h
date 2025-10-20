//____________________________________________________________________________
/*!
\class    genie::DarkSectorDecayer
\brief    Dark Sector decayer module.

          A simple decay simulation...
          ....
          Is a concerete implementation of the EventRecordVisitorI interface.

\author   Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
          University of Sussex

          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  July XX, 2020

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
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


namespace genie {

  class GHepParticle;
  class DarkSectorDecayer : public EventRecordVisitorI {

    using DecayChannel = std::pair<std::vector<int>, double>;
    // first the vector of pdgs, second the decay amplitude

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

    bool ToBeDecayed(const GHepParticle & p) const;
    // comments later
    std::vector<GHepParticle> Decay(const GHepParticle & mother,
                                    const std::vector<int> & pdg_daughters) const;
    // this function will take care of the momentum conservation
    // the output particles cannot be inserted in the event record as they are
    // they need to be translated in space and time, according to the decay amplitude
    int SelectDecayChannel(const std::vector<DecayChannel> & dcs,
                           double total_amplitude) const;
    std::vector<DecayChannel> DarkMediatorDecayChannels(void) const;
    std::vector<DecayChannel> DarkNeutrinoDecayChannels( int mother_pdg ) const;
    void SetSpaceTime(std::vector<GHepParticle> & pp, const GHepParticle & mother,
                      double total_amplitude) const;
    
  private:

    static string ParticleGunKineAsString(const TLorentzVector & vec4 ) ;

    mutable TGenPhaseSpace fPhaseSpaceGenerator;

    double fEps2;
    std::array<double, 4> fMixing2s;
    double fAlpha_D;

    double fDNuMass, fDNuMass2;
    double fDMediatorMass, fDMediatorMass2;

  };

}         // genie namespace
#endif    // _DARK_SECTOR_DECAYER_H_
