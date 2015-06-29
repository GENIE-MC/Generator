//____________________________________________________________________________
/*!

\program gtestMessenger

\brief   Program used for testing / debugging log4cpp

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created June 10, 2004

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"

using namespace genie;

int main(int /*argc*/, char ** /*argv*/)
{
  LOG("Stream-Name", pFATAL)  << "this is a message with priority: FATAL" ;
  LOG("Stream-Name", pALERT)  << "this is a message with priority: ALERT" ;
  LOG("Stream-Name", pCRIT)   << "this is a message with priority: CRIT"  ;
  LOG("Stream-Name", pERROR)  << "this is a message with priority: ERROR" ;
  LOG("Stream-Name", pWARN)   << "this is a message with priority: WARN"  ;
  LOG("Stream-Name", pNOTICE) << "this is a message with priority: NOTICE";
  LOG("Stream-Name", pINFO)   << "this is a message with priority: INFO"  ;
  LOG("Stream-Name", pDEBUG)  << "this is a message with priority: DEBUG" ;

  //-- set a "verbosity level" for the "Stream-Name" message stream
  
  Messenger * msg = Messenger::Instance();
  msg->SetPriorityLevel("Stream-Name", pERROR);

  //-- re-print / some of them must be filtered now
  
  LOG("Stream-Name", pFATAL)  << "this is another message with priority: FATAL" ;
  LOG("Stream-Name", pALERT)  << "this is another message with priority: ALERT" ;
  LOG("Stream-Name", pCRIT)   << "this is another message with priority: CRIT"  ;
  LOG("Stream-Name", pERROR)  << "this is another message with priority: ERROR" ;
  LOG("Stream-Name", pWARN)   << "this is another message with priority: WARN"  ;
  LOG("Stream-Name", pNOTICE) << "this is another message with priority: NOTICE";
  LOG("Stream-Name", pINFO)   << "this is another message with priority: INFO"  ;
  LOG("Stream-Name", pDEBUG)  << "this is another message with priority: DEBUG" ;

  //-- set pririty level to minimum

  msg->SetPriorityLevel("Stream-Name", pDEBUG);
                                
  //-- try different pre-processor macros for printing out messages
  
  LOG_FATAL  ("Stream-Name") << "this is yet another message with priority: FATAL" ;
  LOG_ALERT  ("Stream-Name") << "this is yet another message with priority: ALERT" ;
  LOG_CRIT   ("Stream-Name") << "this is yet another message with priority: CRIT"  ;
  LOG_ERROR  ("Stream-Name") << "this is yet another message with priority: ERROR" ;
  LOG_WARN   ("Stream-Name") << "this is yet another message with priority: WARN"  ;
  LOG_NOTICE ("Stream-Name") << "this is yet another message with priority: NOTICE";
  LOG_INFO   ("Stream-Name") << "this is yet another message with priority: INFO"  ;
  LOG_DEBUG  ("Stream-Name") << "this is yet another message with priority: DEBUG" ;
  
  return 0;
}

