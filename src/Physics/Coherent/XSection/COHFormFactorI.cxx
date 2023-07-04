//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Marco Roda 
         University of Liverpool
         <mroda \at liverpool.ac.uk>
         

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________


#include <limits>

#include "Physics/Coherent/XSection/COHFormFactorI.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

COHFormFactorI::COHFormFactorI( string name ) :
  Algorithm( name ) {
  
}

//____________________________________________________________________________

COHFormFactorI::COHFormFactorI( string name, string config ) :
  Algorithm( name, config ) {
  
}

//____________________________________________________________________________
genie::Range1D_t COHFormFactorI::QRange( int /* pdg */ ) const {
  
  return genie::Range1D_t( 0.,  std::numeric_limits<double>::infinity() ) ;
}
//____________________________________________________________________________
void COHFormFactorI::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHFormFactorI::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
