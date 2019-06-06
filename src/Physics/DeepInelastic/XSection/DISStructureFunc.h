//____________________________________________________________________________
/*!

\class    genie::DISStructureFunc

\brief    A class holding Deep Inelastic Scattering (DIS) Form Factors
          (invariant structure funstions)

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          DISStructureFuncModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 05, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIS_STRUCTURE_FUNCTIONS_H_
#define _DIS_STRUCTURE_FUNCTIONS_H_

#include <iostream>

#include "Physics/DeepInelastic/XSection/DISStructureFuncModelI.h"
#include "Framework/Interaction/Interaction.h"

using std::ostream;

namespace genie {

class DISStructureFunc;
ostream & operator << (ostream & stream, const DISStructureFunc & sf);

class DISStructureFunc {

public:
  DISStructureFunc();
  DISStructureFunc(const DISStructureFunc & form_factors);
  virtual ~DISStructureFunc() { }

  //! Attach an algorithm
  void SetModel  (const DISStructureFuncModelI * model);

  //! Calculate the S/F's for the input interaction using the attached algorithm
  void Calculate (const Interaction * interaction);

  //! Get the computed structure function F1
  double F1 (void) const { return fF1; }

  //! Get the computed structure function F2
  double F2 (void) const { return fF2; }

  //! Get the computed structure function F3
  double F3 (void) const { return fF3; }

  //! Get the computed structure function F4
  double F4 (void) const { return fF4; }

  //! Get the computed structure function F5
  double F5 (void) const { return fF5; }

  //! Get the computed structure function F6
  double F6 (void) const { return fF6; }

  //! Get the attached model
  const DISStructureFuncModelI * Model (void) const {return fModel;}

  void   Reset    (Option_t * opt="");
  void   Copy     (const DISStructureFunc & sf);
  bool   Compare  (const DISStructureFunc & sf) const;
  void   Print    (ostream & stream) const;

  bool               operator == (const DISStructureFunc & sf) const;
  DISStructureFunc & operator =  (const DISStructureFunc & sf);
  friend ostream &   operator << (ostream & stream, const DISStructureFunc & sf);

private:

  double fF1;
  double fF2;
  double fF3;
  double fF4;
  double fF5;
  double fF6;

  const DISStructureFuncModelI * fModel;
};

}       // genie namespace

#endif  // _DIS_STRUCTURE_FUNCTIONS_H_
