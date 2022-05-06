//----------------------------------------------------------------------------
/*!

  Handles the inclusive-type transformation selection
  Assumes pseudounitarity (|U_{e4}|^{2} + |U_{\mu 4}|^{2} = 1) to conserve stats
  and also to be consistent with the fluxes read in from dk2nu
  
  Will have to upscale POT by factor (|U_{e4}|^{2} + |U_{\mu4}^{2}|)^{-2} later!

\namespace  genie::NHL::NHLSelector

\brief      Transformation inclusive-method channel selector

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    December 9th, 2021

\cpright    ??? - TBD

*/
//----------------------------------------------------------------------------

#ifndef _NHL_DECAYSELECTOR_H_
#define _NHL_DECAYSELECTOR_H_

// -- C++ includes
#include <map>

// -- GENIE includes
#include "Framework/Conventions/Constants.h"

#include "Physics/NeutralHeavyLepton/NHLBRFunctions.h"
#include "Physics/NeutralHeavyLepton/NHLDecayMode.h"

namespace genie {
namespace NHL {

    namespace NHLSelector {

      // only need to calculate decay widths once! Store them in this array
      static __attribute__((unused)) double fDecayGammas[] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
      
      // valid channels with widths
      std::map< genie::NHL::NHLDecayMode_t, double > GetValidChannelWidths( const double M, const double Ue42, const double Umu42, const double Ut42, const bool IsMajorana = false );
      // derived
      double GetTotalDecayWidth( std::map< genie::NHL::NHLDecayMode_t, double > gammaMap );
      double CalcCoMLifetime( const double M, const double Ue42, const double Umu42, const double Ut42, const bool IsMajorana = false );
      
      // now choose channels you're interested in seeing
      std::map< genie::NHL::NHLDecayMode_t, double > SetInterestingChannels( std::vector< genie::NHL::NHLDecayMode_t > intChannels, std::map< genie::NHL::NHLDecayMode_t, double > gammaMap );
      
      // transform widths to probabilities (marginalised over all interesting!)
      std::map< genie::NHL::NHLDecayMode_t, double > GetProbabilities( std::map< genie::NHL::NHLDecayMode_t, double > gammaMap );
      
      // make a choice from interesting channels
      // This is the inclusive method - see messages to Xianguo, Dec 9th 2021
      genie::NHL::NHLDecayMode_t SelectChannelInclusive( std::map< genie::NHL::NHLDecayMode_t, double > Pmap, double ranThrow );
      
    } // namespace NHLSelector
    
} // namespace NHL
} // namespace genie

#endif // #ifndef _NHL_DEDCAYSELECTOR_H_
