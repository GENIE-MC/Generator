//____________________________________________________________________________
/*!

\class    genie::COHFormFactorMap

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

#ifndef _COH_FORM_FACTOR_MAP_H_
#define _COH_FORM_FACTOR_MAP_H_

#include <map>

#include "Physics/Coherent/XSection/COHFormFactorI.h"
#include "Physics/Coherent/XSection/DeVriesFormFactor.h"

namespace genie {

class COHFormFactorMap : public COHFormFactorI {

public:

  COHFormFactorMap();
  COHFormFactorMap( string config );
  virtual ~COHFormFactorMap();

  // methods to be implemented from COHFormFactorI
  virtual double ProtonFF ( double Q, int pdg ) const  ;

  virtual double NeutronFF( double Q, int pdg ) const {
    return ProtonFF( Q, pdg ) ;
  }

  virtual bool HasNucleus( int pdg ) const ;

  // methods to implemented from Algorithm
  void Configure (const Registry & config);
  void Configure (string param_set);


 protected:

  COHFormFactorMap( string name, string config );

  virtual void LoadConfig(void);

  std::map<int, const genie::DeVriesFormFactor *> Map() const { return fNuclearFFs; }

private:

  std::map<int, const genie::DeVriesFormFactor *> fNuclearFFs ;
  // the map key is given by the pdg

};

}       // genie namespace
#endif  //  #ifndef _COH_FORM_FACTOR_MAP_H_
