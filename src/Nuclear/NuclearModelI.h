//____________________________________________________________________________
/*!

\class    genie::NuclearModelI

\brief    Pure abstract base class.
          Defines the NuclearModelI interface to be implemented by any physics 
          model describing the distribution of nucleons within a nuclei

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 09, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NUCLEAR_MODEL_I_H_
#define _NUCLEAR_MODEL_I_H_

#include <TVector3.h>

#include "Algorithm/Algorithm.h"
#include "Interaction/Target.h"
#include "Nuclear/NuclearModel.h"

namespace genie {

class NuclearModelI : public Algorithm {

public:
  virtual ~NuclearModelI();

  virtual bool           GenerateNucleon (const Target &) const = 0;
  virtual double         Prob            (double p, double w, const Target &) const = 0;
  virtual NuclearModel_t ModelType       (const Target &) const = 0;

  virtual double         RemovalEnergy   (void)           const;
  virtual double         Momentum        (void)           const;
  virtual TVector3       Momentum3       (void)           const;
  virtual FermiMoverInteractionType_t GetFermiMoverInteractionType(void) const;

protected:
  NuclearModelI();
  NuclearModelI(string name);
  NuclearModelI(string name, string config);

  mutable double   fCurrRemovalEnergy;
  mutable TVector3 fCurrMomentum;
  mutable FermiMoverInteractionType_t fFermiMoverInteractionType;

};

}         // genie namespace
#endif    // _NUCLEAR_MODEL_I_H_

