//____________________________________________________________________________
/*!

\class    genie::AhrensNCELStrangeFF

\brief    Implements the NCELStrangeFormFactorsModelI interface.
\todo     TODO: Cite Ahrens paper

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Feb 23, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _AHRENS_NCEL_STRANGE_FORM_FACTORS_I_H_
#define _AHRENS_NCEL_STRANGE_FORM_FACTORS_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Interaction.h"

#include "Physics/QuasiElastic/XSection/AxialFormFactor.h"
#include "Physics/QuasiElastic/XSection/AxialFormFactorModelI.h"
#include "Physics/QuasiElastic/XSection/NCELStrangeFormFactorsModelI.h"

namespace genie {

class AhrensNCELStrangeFF : public NCELStrangeFormFactorsModelI {

public:
  AhrensNCELStrangeFF();
  AhrensNCELStrangeFF( std::string config );
  virtual ~AhrensNCELStrangeFF();

  // NCELStrangeFormFactorsModelI interface implementation
  inline virtual double F1Vs (const Interaction* /*interaction*/) const
    { return 0.; }
  inline virtual double F2Vs (const Interaction* /*interaction*/) const
    { return 0.; }
  inline virtual double FPs  (const Interaction* /*interaction*/) const
    { return 0.; }

  virtual double FAs(const Interaction* interaction) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  virtual void Configure(const Registry & config);
  virtual void Configure(string config);

  virtual void LoadConfig (void);

protected:

  const AxialFormFactorModelI * fAxFFModel;
  mutable AxialFormFactor fAxFF;

  double fEta;

};

} // genie namespace
#endif // _AHRENS_NCEL_STRANGE_FORM_FACTORS_I_H_
