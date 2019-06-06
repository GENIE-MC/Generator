//____________________________________________________________________________
/*!

\class    genie::COHKinematicsGenerator

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

#ifndef _COH_KINEMATICS_GENERATOR_H_
#define _COH_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

class TF2;

namespace genie {

  class COHKinematicsGenerator : public KineGeneratorWithCache {

  public :
    COHKinematicsGenerator();
    COHKinematicsGenerator(string config);
    ~COHKinematicsGenerator();

    // implement the EventRecordVisitorI interface
    void ProcessEventRecord(GHepRecord * event_rec) const;

    // overload the Algorithm::Configure() methods to load private data
    // members from configuration options
    void Configure(const Registry & config);
    void Configure(string config);

    // methods to load sub-algorithms and config data from the Registry
    void LoadConfig (void);

    // different kinematics calculators for different models
    void   CalculateKin_ReinSehgal(GHepRecord * event_rec) const;
    void   CalculateKin_BergerSehgal(GHepRecord * event_rec) const;
    void   CalculateKin_BergerSehgalFM(GHepRecord * event_rec) const;
    void   CalculateKin_AlvarezRuso(GHepRecord * event_rec) const;
    void SetKinematics(const double E_l, const double theta_l, const double phi_l, 
                       const double theta_pi, const double phi_pi, 
                       const     Interaction* interaction, Kinematics* kinematics) const;
    bool CheckKinematics(const double E_l, const double theta_l, 
                         const double phi_l, const double theta_pi, 
                         const double phi_pi, const Interaction* interaction) const;

    // overload KineGeneratorWithCache method to compute max xsec
    double ComputeMaxXSec (const Interaction * in) const;
    double MaxXSec_ReinSehgal (const Interaction * in) const;
    double MaxXSec_BergerSehgal (const Interaction * in) const;
    double MaxXSec_BergerSehgalFM (const Interaction * in) const;
    double MaxXSec_AlvarezRuso (const Interaction * in) const;

    // overload KineGeneratorWithCache method to get energy
    double Energy         (const Interaction * in) const;

    // TODO: should fEnvelope and fRo be public? They look like they should be private
    mutable TF2 * fEnvelope; ///< 2-D envelope used for importance sampling
    double fRo;              ///< nuclear scale parameter

  private:
    double pionMass(const Interaction* in) const;
    void   throwOnTooManyIterations(unsigned int iters, GHepRecord* evrec) const;

    double fQ2Min;  ///< lower bound of integration for Q^2 in Berger-Sehgal Model
    double fQ2Max;  ///< upper bound of integration for Q^2 in Berger-Sehgal Model
    double fTMax;   ///< upper bound for t = (q - p_pi)^2
  };

}      // genie namespace
#endif // _COH_KINEMATICS_GENERATOR_H_
