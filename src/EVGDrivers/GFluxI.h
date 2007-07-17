//____________________________________________________________________________
/*!

\class    genie::GFluxI

\brief    GENIE Interface for user-defined flux classes

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 25, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _G_FLUX_I_H_
#define _G_FLUX_I_H_

class TLorentzVector;

namespace genie {

class PDGCodeList;

class GFluxI {

public :

  virtual ~GFluxI();

  //-- define the GFluxI interface

  // declare list of neutrinos and maximum energy
  virtual const PDGCodeList &    FluxParticles (void) = 0;
  virtual double                 MaxEnergy     (void) = 0;

  // generate flux neutrino
  virtual bool                   GenerateNext  (void) = 0;
  virtual int                    PdgCode       (void) = 0;
  virtual double                 Weight        (void) = 0;
  virtual const TLorentzVector & Momentum      (void) = 0;
  virtual const TLorentzVector & Position      (void) = 0;

protected:

  GFluxI();
};

}      // genie namespace

#endif // _G_FLUX_I_H_
