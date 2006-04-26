//____________________________________________________________________________
/*!

\class    genie::HadronizationModelI

\brief    Pure abstract base class.
          Defines the HadronizationModelI interface to be implemented by any
          algorithmic class performing hadronization.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 17, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _HADRONIZATION_MODEL_I_H_
#define _HADRONIZATION_MODEL_I_H_

#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "Algorithm/Algorithm.h"
#include "Interaction/Interaction.h"

namespace genie {

class HadronizationModelI : public Algorithm {

public:

  virtual ~HadronizationModelI();

  //-- define HadronizationModelI interface

  virtual void           Initialize   (void)                 const = 0;
  virtual TClonesArray * Hadronize    (const Interaction * ) const = 0;

protected:

  HadronizationModelI();
  HadronizationModelI(string name);
  HadronizationModelI(string name, string config);
};

}         // genie namespace

#endif    // _HADRONIZATION_MODEL_I_H_

