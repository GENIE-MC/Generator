//____________________________________________________________________________
/*!

\class    genie::NuclearModelI

\brief    Pure abstract base class.
          Defines the NuclearModelI interface to be implemented by any physics 
          model describing the distribution of nucleons within a nuclei

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 09, 2004

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
  virtual double         RemovalEnergy   (void)           const;
  virtual double         Momentum        (void)           const;
  virtual TVector3       Momentum3       (void)           const;
  virtual double         Prob            (double p, double w, const Target &) const = 0;
  virtual NuclearModel_t ModelType       (void)           const = 0;

protected:
  NuclearModelI();
  NuclearModelI(string name);
  NuclearModelI(string name, string config);

  mutable double   fCurrRemovalEnergy;
  mutable TVector3 fCurrMomentum;

};

}         // genie namespace
#endif    // _NUCLEAR_MODEL_I_H_

