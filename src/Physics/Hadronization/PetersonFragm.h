//____________________________________________________________________________
/*!

\class    genie::PetersonFragm

\brief    The Peterson fragmentation function.
          Is a concrete implementation of the FragmentationFunctionI interface.

\ref      C.Peterson et al., Phys.Rev.D23, 56 (1981)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  June 15, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PETERSON_FRAGM_H_
#define _PETERSON_FRAGM_H_

#include <TF1.h>

#include "Physics/Hadronization/FragmentationFunctionI.h"

namespace genie {

class PetersonFragm : public FragmentationFunctionI {

public:
  PetersonFragm();
  PetersonFragm(string config);
  ~PetersonFragm();

  //! implement the FragmentationFunctionI interface
  double Value     (double z) const;
  double GenerateZ (void)     const;

  //! methods overloading the Algorithm() interface implementation
  //! to build the fragmentation function from configuration data
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void BuildFunction (void);
  TF1 * fFunc;
};

}      // genie namespace

#endif // _PETERSON_FRAGM_H_
