//____________________________________________________________________________
/*!

\class    genie::DeltaInMediumCorrections

\brief    Transition form factor from Nucleon to Delta
          To be used for the evaluation of the Cross section of 
          COH NC Gamma production


\author  Marco Roda <mroda@liverpool.ac.uk>
         University of Liverpool

\created January 2020

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DELTA_IN_MEDIUM_CORRECTIONS_H_
#define _DELTA_IN_MEDIUM_CORRECTIONS_H_ 

#include "Framework/Algorithm/Algorithm.h"

#include "Physics/NuclearState/FermiMomentumTable.h"


namespace genie {

class DeltaInMediumCorrections : public Algorithm {

public:
  DeltaInMediumCorrections();
  DeltaInMediumCorrections(string config);
  virtual ~DeltaInMediumCorrections();

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

  double Sigma( int pdg ) const ;

  /*
    
    ** Comments about Fermi Momentum and Density ** 

    The corrections to the delta are implemented in the NC COH gamma model as average corrections.
    There are no integrals perforemed at run time. 
    
    Hence, the density is not space dependent and it's an average on the nuclus. 
    It is also the average of the two nucleon densities, that is why there are two methods.

    The model is developed in the context of a Local Fermi Gas model, so tempering with 
    the density means tempering with the density means tempering with the fermi momentum.
    Some degrees of consistency need to be granted.

    Because in these model the nuclear corrections are applied as averages, the FermiMomentum cannot be evaluated as a function of the position, simply because that is not a degree of freedome. 
    So the consistency is achieved in this way:
    1) the fermi momenta for proton and neutrons are taken from the GENIE tables (as if it was a Relativistic Fermi Gas
    2) from the fermi Momenta, we evaluate the average densities. The nucleus density is the average of the two
    3) that average density is what is used to evaluate all the corrections, e.g. the Sigma corrections. 
    This approach should grant consistency within the GENIE framework and with the model itself
    
   */

  double FermiMomentum( int pdg ) const ;
  double FermiMomentum( int nucleus_pdg, int nucleon_pdg ) const ;

  double AverageDensity( int pdg ) const ;
  double AverageDensity( int nucleus_pdg, int nucleon_pdg ) const ;

  double Sigma( int pdg ) const  ;
  
  // AverageDirectPropagator() const ;
  // AverageCrossPropagator() const ;

private:

  void LoadConfig(void);

  const FermiMomentumTable * fKFTable = nullptr ;
  // this object is retrieved with a Pool that keeps ownerships of the tables
  // No need to delete this object 

  double fDeltaV0 ; 
  double fRho0 ;  // this is the nuclear matter density. 
                  // It should be in natural units but it's read from the xml in fm^-3 as that is how it's usually reported

};

}       // genie namespace
#endif  // _DELTA_TRANSITION_FORM_FACTOR_H_  
