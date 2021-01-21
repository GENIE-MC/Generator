//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory 

 Marco Roda <mroda \at liverpool.ac.uk>                                                                    
 University of Liverpool                                                                                   
 
*/
//____________________________________________________________________________

#include "Physics/Common/ParticleLibrary.h"

using namespace genie;


ParticleLibrary::ParticleLibrary( const string & name, const string & config ) :
  Algorithm(name, config) {
  ;
}
//___________________________________________________________________________
ParticleLibrary::~ParticleLibrary() {
  ;
}
//___________________________________________________________________________
void ParticleLibrary::Configure( const Registry & config) {

  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void ParticleLibrary::Configure( string config ) {

  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
