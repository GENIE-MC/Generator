//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelI

\brief    Pure Abstract Base Class. Defines the DISStructureFuncModelI
          interface to be implemented by any algorithmic class computing DIS
          structure functions.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIS_STRUCTURE_FUNCTIONS_MODEL_I_H_
#define _DIS_STRUCTURE_FUNCTIONS_MODEL_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Interaction.h"

namespace genie {

class DISStructureFuncModelI : public Algorithm {

public:
  virtual ~DISStructureFuncModelI();

  //! Calculate the structure functions F1-F6 for the input interaction
  virtual void Calculate (const Interaction * interaction) const = 0;

  //! Get the computed structure function F1
  virtual double F1 (void) const = 0;

  //! Get the computed structure function F2
  virtual double F2 (void) const = 0;

  //! Get the computed structure function F3
  virtual double F3 (void) const = 0;

  //! Get the computed structure function F4
  virtual double F4 (void) const = 0;

  //! Get the computed structure function F5
  virtual double F5 (void) const = 0;

  //! Get the computed structure function F6
  virtual double F6 (void) const = 0;

protected:
  DISStructureFuncModelI();
  DISStructureFuncModelI(string name);
  DISStructureFuncModelI(string name, string config);
};

}         // genie namespace
#endif    // _DIS_STRUCTURE_FUNCTIONS_MODEL_I_H_
