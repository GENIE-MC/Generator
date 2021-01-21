//____________________________________________________________________________
/*!

  \class    genie::PDGParticleLibrary

  \brief    Implementation of the PDG particle library via the ParticleLibrary interface
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

#ifndef _COMMON_PDG_PARTILCE_LIBRARY_H_
#define _COMMON_PDG_PARTILCE_LIBRARY_H_

#include <memory>

#include "Framework/Common/ParticleLibrary.h" 

#include "TDatabasePDG.h"


namespace genie {

  class PDGParticleLibrary : public ParticleLibrary {

  public :
    
    PDGParticleLibrary() ;
    PDGParticleLibrary( const string & config ) ;
    PDGParticleLibrary( const PDGParticleLibrary & ) = delete ;
    PDGParticleLibrary( PDGParticleLibrary && ) = delete ;

    virtual ~PDGParticleLibrary() ;

    // Interfaces defined by PDGParticleLibrary
    virtual TParticlePDG * Find(int pdgc) const override ; ///< The TParticlPDG returned here is owned by the Library
                                                   
  protected:
    
    virtual void LoadConfig( void ) override ; 
  
    TDatabasePDG & Database() { return *fDatabase ; }

  private :

    std::unique_ptr<TDatabasePDG> fDatabase ;

    std::string PDGTableFile() const ;

  };

}      // genie namespace

#endif // _COMMON_PDG_PARTILCE_LIBRARY_H_
