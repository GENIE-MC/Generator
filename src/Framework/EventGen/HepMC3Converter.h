//____________________________________________________________________________
/*!

\class    genie::HepMC3Converter

\brief    Converts genie::EventRecord objects into HepMC3::GenEvent objects
          following the NuHepMC standard (https://github.com/NuHepMC/Spec)

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  January 3, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HEPMC3_CONVERTER_H_
#define _HEPMC3_CONVERTER_H_

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_HEPMC3_INTERFACE_ENABLED__
// Forward-declare needed HepMC3 classes here
namespace HepMC3 {
  class GenEvent;
  class GenRunInfo;
}

namespace genie {

class EventRecord;
class GMCJDriver;

class HepMC3Converter {

public:

  HepMC3Converter(void);

  std::shared_ptr< HepMC3::GenEvent > ConvertToHepMC3(
    const genie::EventRecord& gevrec );

  std::shared_ptr< genie::EventRecord > RetrieveGHEP(
    const HepMC3::GenEvent& evt );

  void AttachGMCJDriver( const genie::GMCJDriver* mc_driver );

protected:

  int GetNuHepMCParticleStatus( const genie::GHepParticle* gpart,
    const genie::EventRecord& gevrec ) const;

  int GetNuHepMCProcessID( const genie::Interaction& inter ) const;

  genie::GHepStatus_t GetGHepParticleStatus( int nuhepmc_status ) const;

  void StoreInteraction( const genie::Interaction& inter,
    HepMC3::GenEvent& evt );

  genie::Interaction* RetrieveInteraction( const HepMC3::GenEvent& evt );

  void PrepareRunInfo( const genie::EventRecord* gevrec );

  void PrepareMCDriverEventInfo( HepMC3::GenEvent& evt );

  std::shared_ptr< HepMC3::GenRunInfo > fRunInfo;

  const genie::GMCJDriver* fMCDriver = nullptr;

};

} // genie namespace

#endif // __GENIE_HEPMC3_INTERFACE_ENABLED__
#endif // _HEPMC3_CONVERTER_H_
