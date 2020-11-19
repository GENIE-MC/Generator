//____________________________________________________________________________
/*!

\class    genie::COHProtonFormFactorInterpolation

\brief    for COH Production Form Factor Model
          The class is develope specifically for the NC COH Gamma
          But in principle these Form Factors could be reused.
          It extends the functionality of COHFormFactorMap providing
          interpolation for the nuclei that are missing from the
          DeVriesFormFactor paper.
          We only interpolate the proton component in this class
          The neutron component will be the proton of the same nuclei rescaled

\author   Marco Roda <mroda@liverpool.ac.uk>
          University of Liverpool

\created  November 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _COH_PROTON_FORM_FACTOR_INTERPOLATION_H_
#define _COH_PROTON_FORM_FACTOR_INTERPOLATION_H_

#include <map>
#include <functional>

#include "Physics/Coherent/XSection/COHFormFactorI.h"

namespace genie {

class COHProtonFormFactorInterpolation : public COHFormFactorI {

public:

  COHProtonFormFactorInterpolation();
  COHProtonFormFactorInterpolation( string config );
  virtual ~COHProtonFormFactorInterpolation();

  // methods to be implemented from COHFormFactorI
  virtual double ProtonFF ( double Q, int pdg ) const override ;
  virtual double NeutronFF( double Q, int pdg ) const override ;

  virtual bool HasNucleus( int pdg ) const override ;

  virtual genie::Range1D_t QRange( int pdg ) const override ;

 protected:

  virtual void LoadConfig(void) override ;

  std::vector<int> Neighbours( int pdg ) const ; 
  
 private:
  
  std::map<int, std::map<int, int>> fArchive ; 
  // the archive is the list of available nuclei in the base Form Factors that you interpolate
  // it is organised with the first key being Z, and the second begin N
  // The contained object is the pdg 

  const genie::COHFormFactorI * fBaseFF ; 
  
  bool fAllowExtrapolation ;
  
};

}       // genie namespace
#endif  //  #ifndef _COH_PROTON_FORM_FACTOR_INTERPOLATION_H_ 
