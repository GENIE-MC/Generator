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

 Important revisions after version 2.0.0 :
 @ Mar 18, 2016 - JJ (SD)
   Added option for GenerateNucleon() to be called with a target and a radius
   as the arguments. Currently used by LocalFGM. Calls 
   GenerateNucleon() with the radius set to 0 for all other NuclearModelI
   implementations.

*/
//____________________________________________________________________________

#ifndef _NUCLEAR_MODEL_I_H_
#define _NUCLEAR_MODEL_I_H_

#include <string>

#include <TVector3.h>

#include "Types/NuclearModel.h"
#include "Algorithm/Algorithm.h"
#include "Interaction/Target.h"

namespace genie {

class NuclearModelI : public Algorithm {

public:
  virtual ~NuclearModelI() {};

  virtual bool           GenerateNucleon (const Target &) const = 0;
  virtual double         Prob            (double p, double w, const Target &) const = 0;
  virtual NuclearModel_t ModelType       (const Target &) const = 0;

  virtual double         RemovalEnergy   (void)           const
  {
    return fCurrRemovalEnergy;
  }
  virtual double         Momentum        (void)           const
  {
    return fCurrMomentum.Mag();
  };
  virtual TVector3       Momentum3       (void)           const
  {
    return fCurrMomentum;
  };
  virtual FermiMoverInteractionType_t GetFermiMoverInteractionType(void) const
  {
    return fFermiMoverInteractionType;
  };

  virtual bool GenerateNucleon(const Target & tgt, 
			       double hitNucleonRadius) const
  {
    return GenerateNucleon(tgt);
  }
  virtual double Prob(double p, double w, const Target & tgt,
	      double hitNucleonRadius) const
  {
    return Prob(p,w,tgt);
  }

protected:
  NuclearModelI()
    : Algorithm(), fFermiMoverInteractionType(kFermiMoveDefault)
    {};
  NuclearModelI(std::string name)
    : Algorithm(name), fFermiMoverInteractionType(kFermiMoveDefault)
    {};
  NuclearModelI(std::string name, std::string config)
    : Algorithm(name, config), fFermiMoverInteractionType(kFermiMoveDefault)
    {};

  mutable double   fCurrRemovalEnergy;
  mutable TVector3 fCurrMomentum;
  mutable FermiMoverInteractionType_t fFermiMoverInteractionType;

};

}         // genie namespace
#endif    // _NUCLEAR_MODEL_I_H_

