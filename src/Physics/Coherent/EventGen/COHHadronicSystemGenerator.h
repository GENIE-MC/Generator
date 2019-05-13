//____________________________________________________________________________
/*!

\class    genie::COHHadronicSystemGenerator

\brief    Generates the f/s hadronic system in v COH pi production interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COH_HADRONIC_SYSTEM_GENERATOR_H_
#define _COH_HADRONIC_SYSTEM_GENERATOR_H_

#include "Physics/Common/HadronicSystemGenerator.h"

namespace genie {

  class XclsTag;

  class COHHadronicSystemGenerator : public HadronicSystemGenerator {

  public :
    COHHadronicSystemGenerator();
    COHHadronicSystemGenerator(string config);
    ~COHHadronicSystemGenerator();

    // implement the EventRecordVisitorI interface
    void ProcessEventRecord(GHepRecord * event_rec) const;
    void CalculateHadronicSystem_ReinSehgal(GHepRecord * event_rec) const;
    void CalculateHadronicSystem_BergerSehgal(GHepRecord * event_rec) const;
    void CalculateHadronicSystem_BergerSehgalFM(GHepRecord * event_rec) const;
    void CalculateHadronicSystem_AlvarezRuso(GHepRecord * event_rec) const;

  private:
    int getPionPDGCodeFromXclTag(const XclsTag& xcls_tag) const;
  };

}      // genie namespace
#endif // _COH_HADRONIC_SYSTEM_GENERATOR_H_

