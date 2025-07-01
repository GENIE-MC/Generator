//----------------------------------------------------------------------------
/*!

  Handles the inclusive-type transformation selection

\namespace  genie::hnl::selector

\brief      Transformation inclusive-method channel selector

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    December 9th, 2021

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org

*/
//----------------------------------------------------------------------------

#ifndef _HNL_DECAYSELECTOR_H_
#define _HNL_DECAYSELECTOR_H_

// -- C++ includes
#include <map>

// -- GENIE includes
#include "Framework/Conventions/Constants.h"

#include "Physics/BeamHNL/HNLBRCalculator.h"
#include "Physics/BeamHNL/HNLDecayMode.h"

namespace genie {
namespace hnl {

    namespace selector {

      // only need to calculate decay widths once! Store them in this array
      static __attribute__((unused)) double fDecayGammas[] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
      
      // valid channels with widths
      std::map< genie::hnl::HNLDecayMode_t, double > GetValidChannelWidths( const double M, const double Ue42, const double Umu42, const double Ut42, const bool IsMajorana = false );
      // derived
      double GetTotalDecayWidth( std::map< genie::hnl::HNLDecayMode_t, double > gammaMap );
      double CalcCoMLifetime( const double M, const double Ue42, const double Umu42, const double Ut42, const bool IsMajorana = false );
      
      // now choose channels you're interested in seeing
      std::map< genie::hnl::HNLDecayMode_t, double > SetInterestingChannels( std::vector< genie::hnl::HNLDecayMode_t > intChannels, std::map< genie::hnl::HNLDecayMode_t, double > gammaMap );
      
      // transform widths to probabilities (marginalised over all interesting!)
      std::map< genie::hnl::HNLDecayMode_t, double > GetProbabilities( std::map< genie::hnl::HNLDecayMode_t, double > gammaMap );
      
      // make a choice from interesting channels
      // This is the inclusive method - see messages to Xianguo, Dec 9th 2021
      genie::hnl::HNLDecayMode_t SelectChannelInclusive( std::map< genie::hnl::HNLDecayMode_t, double > Pmap, double ranThrow );
      
    } // namespace selector
    
} // namespace hnl
} // namespace genie

#endif // #ifndef _HNL_DEDCAYSELECTOR_H_
