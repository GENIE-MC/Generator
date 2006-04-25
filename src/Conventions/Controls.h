//____________________________________________________________________________
/*!

\namespace genie::controls

\brief     Misc GENIE control constants

\author    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
           CCLRC, Rutherford Appleton Laboratory

\created   May 03, 2004

*/
//____________________________________________________________________________

#ifndef _CONTROLS_H_
#define _CONTROLS_H_

namespace genie {
namespace controls {

// maximum allowed number of iterations in rejection MC method
// before selecting a valid number
static const unsigned int kRjMaxIterations = 1000;

// maximum allowed depth when GENIE is running in recursive mode
static const unsigned int kRecursiveModeMaxDepth = 100;

// maximum allowed number of EVGThreadExceptions that is allowed
// to be caught by EventGenerator at a single event generation thread
static const unsigned int kMaxEVGThreadExceptions = 350;

// Default random number generator seed number. It can be overriden
// setting the $GSEED env. var. or by using RandomGen::SetSeed(int)
static const unsigned int kDefaultRandSeed = 65539;

//----- Misc hard limits, cuts
static const double kMinQ2Limit   = 1e-5;  // GeV^2
static const int kMaxMultiplicity = 35;    // for KNO hadronization model


} // namespace controls
} // namespace genie

#endif // _CONTROLS_H_


