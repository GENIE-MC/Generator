//____________________________________________________________________________
/*!

\class    genie::NHLDecayMode

\brief    Enumeration of NHL decay modes.

\author

\created  November 10, 2011

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NHL_DECAY_MODE_H_
#define _NHL_DECAY_MODE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

 typedef enum ENHLDecayMode {

   kNHLDcyNull = 0,

   //
   // Add enumerations for each decay channel
   // Start eacg enum with kNHLDcy...
   // See $GENIE/src/Physics/NucleonDecay/NucleonDecayMode.h for a similar code
   //
   //
   // ...
   // ...
   // ...
   kNHLDcyTEST = 1


 } NHLDecayMode_t;

}
#endif
