//____________________________________________________________________________
/*!

\namespace genie::constants::intranuke

\brief     Data used by INTRANUKE cascade MC for intranuclear rescattering.

\author    name <e-mail>

\created   Month xx, yyyy

*/
//____________________________________________________________________________

#ifndef _INTRANUKE_CONSTANTS_H_
#define _INTRANUKE_CONSTANTS_H_

namespace genie     {
namespace constants {
namespace intranuke {

// Cummulative interaction probabilities for pi+Fe in 50 MeV bins
// Data from NeuGEN's Intranuke

static const int kPNDataPoints = 12;

static const double kPElastic[kPNDataPoints] = {
       .97,.94,.93,.92,.91,.90,.90,.90,.90,.90,.90,.90
};

static const double kPInelastic[kPNDataPoints] = {
       .51,.53,.55,.56,.56,.57,.57,.57,.57,.57,.57,.57
};

static const double kPAbsorption[kPNDataPoints] = {
       .13,.18,.23,.23,.18,.13,.10,.08,.07,.06,.05,.05
};

} // intranuke namespace
} // constants namespace
} // genie     namespace

#endif // _INTRANUKE_CONSTANTS_H_
