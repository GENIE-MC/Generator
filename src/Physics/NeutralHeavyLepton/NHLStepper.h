//----------------------------------------------------------------------------
/*!

  Steps NHL forward by time-step and determines if it's decayed yet
  If inclusive: just pick position of vertex
  If exclusive: also query which channel ran out first

\namespace  genie::NHL::NHLSelector

\brief      Vertex position calculator

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    December 10th, 2021

\cpright    ??? - TBD

*/
//----------------------------------------------------------------------------

#ifndef _NHL_STEPPER_H_
#define _NHL_STEPPER_H_

// -- GENIE includes

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
  
#include "Physics/NeutralHeavyLepton/SimpleNHL.h"

namespace genie {
namespace NHL {

    class SimpleNHL;

    namespace NHLSelector {
	
      extern std::map< genie::NHL::NHLDecayMode_t, double > fCounterMap;
      const double fdt = 0.01 * genie::units::nanosecond; // in lab frame
      
      // inclusive method stepper
      std::vector< double > PropagateTilDecay( genie::NHL::SimpleNHL sh ); // returns decay 4-vertex
      
      // exclusive method stepper
      void PropagateAndSelectChannel( genie::NHL::SimpleNHL sh );
      
      // empty the vector of counters
      void emptyCounters( );
      
      // populate the vector of counters
      void fillCounters( std::map< genie::NHL::NHLDecayMode_t, double > gammaMap );

    } // namespace NHLSelector
    
} // namespace NHL
} // namespace genie

#endif // #ifndef _NHL_STEPPER_H_
