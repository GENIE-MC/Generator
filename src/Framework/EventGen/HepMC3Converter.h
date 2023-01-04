//____________________________________________________________________________
/*!

\class    genie::HepMC3Converter

\brief    Converts genie::EventRecord objects into HepMC3::GenEvent objects
          following the NuHepMC standard (https://github.com/NuHepMC/Spec)

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  January 3, 2023

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
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
  class WriterAscii;
}

namespace genie {

class EventRecord;

class HepMC3Converter {

public:

  HepMC3Converter(void);

  std::shared_ptr< HepMC3::GenEvent > ConvertToHepMC3(
    const genie::EventRecord& gevrec );

  // For testing purposes
  // TODO: Remove and create an ntuple writer instead
  void WriteEvent( const HepMC3::GenEvent& evt ) const;

protected:

  int GetHepMC3ParticleStatus( const genie::GHepParticle* gpart,
    const genie::EventRecord& gevrec ) const;

  void PrepareRunInfo();

  std::shared_ptr< HepMC3::GenRunInfo > fRunInfo;
  std::shared_ptr< HepMC3::WriterAscii > fWriter;
};

} // genie namespace

#endif // __GENIE_HEPMC3_INTERFACE_ENABLED__
#endif // _HEPMC3_CONVERTER_H_
