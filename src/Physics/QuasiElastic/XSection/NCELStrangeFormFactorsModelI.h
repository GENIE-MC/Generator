//____________________________________________________________________________
/*!

\class    genie::NCELStrangeFormFactorsModelI

\brief    Pure abstract base class. Defines the NCELStrangeFormFactorsModelI
          interface to be implemented by any algorithmic class computing
          strange form factors for neutral current elastic scattering

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Feb 23, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NCEL_STRANGE_FORM_FACTORS_MODEL_I_H_
#define _NCEL_STRANGE_FORM_FACTORS_MODEL_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Interaction.h"

namespace genie {

class NCELStrangeFormFactorsModelI : public Algorithm {

public:
  virtual ~NCELStrangeFormFactorsModelI();

  //! Compute the strange form factor F1Vs for the input interaction
  virtual double F1Vs(const Interaction* interaction) const = 0;

  //! Compute the strange form factor F2Vs for the input interaction
  virtual double F2Vs(const Interaction* interaction) const = 0;

  //! Compute the strange form factor FAs for the input interaction
  virtual double FAs(const Interaction* interaction) const = 0;

  //! Compute the strange form factor FPs for the input interaction
  virtual double FPs(const Interaction* interaction) const = 0;

  //! \name Second-class currents
  //! @{
  //! Compute the strange form factor F3Vs for the input interaction
  inline virtual double F3Vs(const Interaction* /*interaction*/)
    const { return 0.; }

  //! Compute the strange form factor F3As for the input interaction
  inline virtual double F3As(const Interaction* /*interaction*/)
    const { return 0.; }
  //! @}

protected:
  NCELStrangeFormFactorsModelI();
  NCELStrangeFormFactorsModelI( std::string name );
  NCELStrangeFormFactorsModelI( std::string name, std::string config );
};

} // genie namespace
#endif // _NCEL_STRANGE_FORM_FACTORS_MODEL_I_H_
