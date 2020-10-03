//____________________________________________________________________________
/*!

\class    genie::COHGammaKinematicsGenerator

\brief    Generates values for the kinematic variables describing coherent 
          neutrino-nucleus pion production events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COH_GAMMA_KINEMATICS_GENERATOR_H_
#define _COH_GAMMA_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

#include "Physics/Coherent/XSection/COHGammaIntegrationLimits.h" 

namespace genie {

  class COHGammaKinematicsGenerator : public KineGeneratorWithCache {

  public :
    COHGammaKinematicsGenerator();
    COHGammaKinematicsGenerator(string config);
    ~COHGammaKinematicsGenerator();

    // implement the EventRecordVisitorI interface
    void ProcessEventRecord(GHepRecord * event_rec) const override ;

    // overload the Algorithm::Configure() methods to load private data
    // members from configuration options
    void Configure(const Registry & config) override ;
    void Configure(string config) override ;

    // methods to load sub-algorithms and config data from the Registry
    void LoadConfig (void);

    // overload KineGeneratorWithCache method to compute max xsec
    double ComputeMaxXSec (const Interaction * in) const override ;

    // overload KineGeneratorWithCache method to get energy
    // This can probably go
    double Energy         (const Interaction * in) const override ;


  protected:
    void   throwOnTooManyIterations(unsigned int iters, GHepRecord* evrec) const;

  private:

    const COHGammaIntegrationLimits * fGammaLimits ; 
    std::array<double, 4> fMinimInitialRatio ;


  };
  
}      // genie namespace
#endif // _COH_GAMMA_KINEMATICS_GENERATOR_H_
