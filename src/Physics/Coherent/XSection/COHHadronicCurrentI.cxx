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


#include "Physics/Coherent/XSection/COHHadronicCurrentI.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

COHHadronicCurrentI::COHHadronicCurrentI( string name ) :
  Algorithm( name ) {
  
}

//____________________________________________________________________________

COHHadronicCurrentI::COHHadronicCurrentI( string name, string config ) :
  Algorithm( name, config ) {
  
}
