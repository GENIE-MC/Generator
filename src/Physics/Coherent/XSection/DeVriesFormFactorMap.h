//____________________________________________________________________________
/*!

\class    genie::DeVriesFormFactorMap

\brief    for COH Production Form Factor Model
          The class is develope specifically for the NC COH Gamma
          But in principle these Form Factors could be reused.

\author   Marco Roda <mroda@liverpool.ac.uk>
          University of Liverpool

\created  November 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DEVRIES_FORM_FACTOR_MAP_H_
#define _DEVRIES_FORM_FACTOR_MAP_H_

#include <map>

#include "Physics/Coherent/XSection/COHFormFactorI.h"
#include "Physics/Coherent/XSection/DeVriesFormFactor.h"

namespace genie {

class DeVriesFormFactorMap : public COHFormFactorI {

public:

  DeVriesFormFactorMap();
  DeVriesFormFactorMap( string config );
  virtual ~DeVriesFormFactorMap();

  // methods to be implemented from COHFormFactorI
  virtual double ProtonFF ( double Q, int pdg ) const override ;

  virtual double NeutronFF( double Q, int pdg ) const override {
    return ProtonFF( Q, pdg ) ;
  }

  virtual bool HasNucleus( int pdg ) const override ;

  virtual genie::Range1D_t QRange( int pdg ) const override ;

  // methods to implemented from Algorithm
  void Configure (const Registry & config) override ;
  void Configure (string param_set) override;

 protected:

  DeVriesFormFactorMap( string name, string config );

  virtual void LoadConfig(void);

  const std::map<int, const genie::DeVriesFormFactor *> & Map() const noexcept { return fNuclearFFs; }

private:

  std::map<int, const genie::DeVriesFormFactor *> fNuclearFFs ;
  // the map key is given by the pdg

};

}       // genie namespace
#endif  //  #ifndef _DEVRIES_FORM_FACTOR_MAP_H_
