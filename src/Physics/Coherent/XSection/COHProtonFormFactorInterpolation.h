//____________________________________________________________________________
/*!

\class    genie::COHProtonFormFactorInterpolation

\brief    for COH Production Form Factor Model
          The class is develope specifically for the NC COH Gamma
          But in principle these Form Factors could be reused.
          It extends the functionality of a generic COHFormFactorI providing
          interpolation for the nuclei that are missing from the 
	  actual implementationthe.
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

  using neutron_map = std::map<int, int> ; 
  // first number of neutron 
  // second pdg

public:

  COHProtonFormFactorInterpolation();
  COHProtonFormFactorInterpolation( string config );
  virtual ~COHProtonFormFactorInterpolation();

  // methods to be implemented from COHFormFactorI
  virtual double ProtonFF ( double Q, int pdg ) const override ;
  virtual double NeutronFF( double Q, int pdg ) const override ;

  virtual bool HasNucleus( int pdg ) const override ;

  virtual genie::Range1D_t QRange( int pdg ) const override ;
  // In this case, it will return the smallest interval that contains all the 
  // ranges from the neighbours used in the interpolation

 protected:

  virtual void LoadConfig(void) ;

  std::vector<int> Neighbours( int pdg ) const ; 

  int ClosestIsotope( const neutron_map & map , int n ) const ;
  // return the closest isotope with than n

  double Interpolate( const vector<int> & zs, const vector<double> & ffs, 
		      int final_z ) const  ;
  // The interpolation interface is done with a generic vector because we think about 
  // more than linear interpolations
    
 private:

  std::map<int, neutron_map> fArchive ; 
  // the archive is the list of available nuclei in the base Form Factors that you interpolate
  // it is organised with the first key being Z, and the second begin N
  // The contained object is the pdg 

  const genie::COHFormFactorI * fBaseFF ; 
  
  bool fAllowExtrapolation ;
  
};

}       // genie namespace
#endif  //  #ifndef _COH_PROTON_FORM_FACTOR_INTERPOLATION_H_ 
