//____________________________________________________________________________
/*!

\class    genie::DeVriesFormFactorInterpolation

\brief    for COH Production Form Factor Model
          The class is develope specifically for the NC COH Gamma
          But in principle these Form Factors could be reused.
          It extends the functionality of COHFormFactorMap providing
          interpolation for the nuclei that are missing from the
          DeVriesFormFactor paper.

\author   Marco Roda <mroda@liverpool.ac.uk>
          University of Liverpool

\created  July 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DEVRIES_FORM_FACTOR_INTERPOLATION_H_
#define _DEVRIES_FORM_FACTOR_INTERPOLATION_H_

#include <map>
#include <functional>

#include "Physics/Coherent/XSection/DeVriesFormFactorMap.h"

namespace genie {

class DeVriesFormFactorInterpolation : public DeVriesFormFactorMap {

public:

  DeVriesFormFactorInterpolation();
  DeVriesFormFactorInterpolation( string config );
  virtual ~DeVriesFormFactorInterpolation();

  // methods to be implemented from COHFormFactorI
  virtual double ProtonFF ( double Q, int pdg ) const override ;
  virtual double NeutronFF( double Q, int pdg ) const override ;

  virtual bool HasNucleus( int pdg ) const override ;

  virtual genie::Range1D_t QRange( int pdg ) const override ;

 protected:

  virtual void LoadConfig(void) override ;

  virtual const genie::FourierBesselFFCalculator & InterpolateProtons( int pdg ) const ;
  virtual const genie::FourierBesselFFCalculator & InterpolateNeutrons( int pdg ) const ;

private:

  pair<int, int> NearbyNuclei( int pdg ) const ;
  double RadiusInterpolation( int pdg, const pair<int, int> & neighbours ) const ;

  genie::FourierBesselFFCalculator LinearInterpolation( int pdg,
                                                        const std::function<int(int)> & ) const ;

  mutable std::map<int, genie::FourierBesselFFCalculator> fInterProtons ;
  mutable std::map<int, genie::FourierBesselFFCalculator> fInterNeutrons ;
  // the map key is given by the pdg

  bool fAllowExtrapolation ;
  
};

}       // genie namespace
#endif  //  #ifndef _DEVRIES_FORM_FACTOR_MAP_H_
