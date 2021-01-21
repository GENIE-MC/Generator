//____________________________________________________________________________
/*!

  \class    genie::ParticleLibrary

  \brief    Interface for a particle library interface
  This allows co-existence of multiple PDG values inside our algorithm

  \author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
            University of Liverpool & STFC Rutherford Appleton Laboratory

            Marco Roda <mroda \at liverpool.ac.uk>
            University of Liverpool

  \created  January 21, 2021

  \cpright  Copyright (c) 2003-2020, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _COMMON_PARTILCE_LIBRARY_H_
#define _COMMON_PARTILCE_LIBRARY_H_

#include "TParticlePDG.h"

#include "Framework/Algorithm/Algorithm.h" 

namespace genie {

  class ParticleLibrary : public Algorithm {

  public :

    ParticleLibrary() = delete ;
    ParticleLibrary( const ParticleLibrary & ) = delete ;
    ParticleLibrary( ParticleLibrary && ) = delete ;

    virtual ~ParticleLibrary() ;

    // Algorithm interfaces facility

    virtual void Configure (const Registry & config) override ;
    virtual void Configure (string config) override ;

    // Interfaces defined by ParticleLibrary
    virtual TParticlePDG * Find(int pdgc) const = 0 ; ///< The TParticlPDG returned here is owned by the Library
                                                   
  protected:
    
    ParticleLibrary( const string & name, const string & config ) ;
    
    virtual void LoadConfig( void ) = 0 ; 
  
  private :

    
  };

}      // genie namespace

#endif // _COMMON_PARTILCE_LIBRARY_H_
