//----------------------------------------------------------------------------
/*!

  Handles the inclusive-type transformation selection
  Assumes pseudounitarity (|U_{e4}|^{2} + |U_{\mu 4}|^{2} = 1) to conserve stats
  and also to be consistent with the fluxes read in from dk2nu
  
  Will have to upscale POT by factor (|U_{e4}|^{2} + |U_{\mu4}^{2}|)^{-2} later!

\namespace  genie::HNL::HNLSelector

\brief      Transformation inclusive-method channel selector

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    December 9th, 2021

\cpright    ??? - TBD

*/
//----------------------------------------------------------------------------

#ifndef _HNL_DECAYSELECTOR_H_
#define _HNL_DECAYSELECTOR_H_

// -- C++ includes
#include <map>

// -- GENIE includes
#include "Framework/Conventions/Constants.h"

#include "Physics/BeamHNL/HNLBRFunctions.h"
#include "Physics/BeamHNL/HNLDecayMode.h"

namespace genie {
namespace HNL {

    namespace HNLSelector {

      // only need to calculate decay widths once! Store them in this array
      static __attribute__((unused)) double fDecayGammas[] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
      
      // valid channels with widths
      std::map< genie::HNL::HNLDecayMode_t, double > GetValidChannelWidths( const double M, const double Ue42, const double Umu42, const double Ut42, const bool IsMajorana = false );
      // derived
      double GetTotalDecayWidth( std::map< genie::HNL::HNLDecayMode_t, double > gammaMap );
      double CalcCoMLifetime( const double M, const double Ue42, const double Umu42, const double Ut42, const bool IsMajorana = false );
      
      // now choose channels you're interested in seeing
      std::map< genie::HNL::HNLDecayMode_t, double > SetInterestingChannels( std::vector< genie::HNL::HNLDecayMode_t > intChannels, std::map< genie::HNL::HNLDecayMode_t, double > gammaMap );
      
      // transform widths to probabilities (marginalised over all interesting!)
      std::map< genie::HNL::HNLDecayMode_t, double > GetProbabilities( std::map< genie::HNL::HNLDecayMode_t, double > gammaMap );
      
      // make a choice from interesting channels
      // This is the inclusive method - see messages to Xianguo, Dec 9th 2021
      genie::HNL::HNLDecayMode_t SelectChannelInclusive( std::map< genie::HNL::HNLDecayMode_t, double > Pmap, double ranThrow );
      
    } // namespace HNLSelector
    
} // namespace HNL
} // namespace genie

#endif // #ifndef _HNL_DEDCAYSELECTOR_H_
