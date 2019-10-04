//____________________________________________________________________________
/*!

\class    genie::GFluxI

\brief    GENIE Interface for user-defined flux classes

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 25, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

 @ Feb 22, 2011 - JD
   Changed the interface by adding the GFluxI::Clear, GFluxI::Index and 
   GFluxI::GenerateWeighted methods needed so that can be used with the new 
   pre-generation of flux interaction probabilities functionality added to
   GMCJDriver. 

*/
//____________________________________________________________________________

#ifndef _G_FLUX_I_H_
#define _G_FLUX_I_H_

#include <TObject.h>

class TLorentzVector;

namespace genie {

class PDGCodeList;

class GFluxI {

public :
  virtual ~GFluxI();

  //
  // define the GFluxI interface:
  //
  virtual const PDGCodeList &    FluxParticles (void) = 0; ///< declare list of flux neutrinos that can be generated (for init. purposes)
  virtual double                 MaxEnergy     (void) = 0; ///< declare the max flux neutrino energy that can be generated (for init. purposes)
  virtual bool                   GenerateNext  (void) = 0; ///< generate the next flux neutrino (return false in err)
  virtual int                    PdgCode       (void) = 0; ///< returns the flux neutrino pdg code
  virtual double                 Weight        (void) = 0; ///< returns the flux neutrino weight (if any)
  virtual const TLorentzVector & Momentum      (void) = 0; ///< returns the flux neutrino 4-momentum 
  virtual const TLorentzVector & Position      (void) = 0; ///< returns the flux neutrino 4-position (note: expect SI rather than physical units)
  virtual bool                   End           (void) = 0; ///< true if no more flux nu's can be thrown (eg reaching end of beam sim ntuples)
  virtual long int               Index         (void) = 0; ///< returns corresponding index for current flux neutrino (e.g. for a flux ntuple returns the current entry number)
  virtual void                   Clear            (Option_t * opt   ) = 0; ///< reset state variables based on opt
  virtual void                   GenerateWeighted (bool gen_weighted) = 0; ///< set whether to generate weighted or unweighted neutrinos

protected:
  GFluxI();
};

}      // genie namespace
#endif // _G_FLUX_I_H_
