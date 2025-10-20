//____________________________________________________________________________
/*!

\class   genie::hnl::ChannelCalculatorI

\brief   Pure abstract base class. Defines the ChannelCalculatorI interface
         to be implemented by BRCalculator Algorithm for calculating HNL production
         and decay rates.

\author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
         University of Oxford

	 based off AxialFormFactorModelI by
	 Aarom Meyer <asmeyer2012 \at uchicago.edu>

	 Costas Andreopoulos <c.andreopoulos \at cern.ch>
	 University of Liverpool

\created November 17, 2022

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HNL_CHANNEL_CALCULATOR_I_H_
#define _HNL_CHANNEL_CALCULATOR_I_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Physics/BeamHNL/HNLDecayMode.h"
#include "Physics/BeamHNL/HNLProductionMode.h"

namespace genie {

  class Registry;

  namespace hnl {

    class ChannelCalculatorI : public Algorithm {

    public:

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      virtual void Configure(const Registry & config) = 0;
      virtual void Configure(string config) = 0;

      // return the kinematic scaling for a production channel
      virtual double KinematicScaling( HNLProd_t hnlprod ) const = 0;

      // return the integrated decay width for a decay channel
      virtual double DecayWidth( HNLDecayMode_t hnldm ) const = 0;

    protected:

      ChannelCalculatorI();
      ChannelCalculatorI(string name);
      ChannelCalculatorI(string name, string config);

    }; // class ChannelCalculatorI
  } // namespace hnl
} // namespace genie

#endif // #ifndef _HNL_CHANNEL_CALCULATOR_I_H_
