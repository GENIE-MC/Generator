//____________________________________________________________________________
/*!

\class    genie::HadronizationModelI

\brief    Pure abstract base class.
          Defines the HadronizationModelI interface to be implemented by any
          algorithmic class performing hadronization.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  August 17, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HADRONIZATION_MODEL_I_H_
#define _HADRONIZATION_MODEL_I_H_

#include "Framework/Algorithm/Algorithm.h"

class TClonesArray;
class TH1D;

namespace genie {

class Interaction;
class PDGCodeList;

class HadronizationModelI : public Algorithm {

public:

  virtual ~HadronizationModelI();

  //! define the HadronizationModelI interface

  virtual void           Initialize       (void)                                 const = 0;
  virtual TClonesArray * Hadronize        (const Interaction * )                 const = 0;
  virtual double         Weight           (void)                                 const = 0;
  virtual PDGCodeList *  SelectParticles  (const Interaction*)                   const = 0;
  virtual TH1D *         MultiplicityProb (const Interaction*, Option_t* opt="") const = 0;

protected:

  HadronizationModelI();
  HadronizationModelI(string name);
  HadronizationModelI(string name, string config);
};

}         // genie namespace

#endif    // _HADRONIZATION_MODEL_I_H_

