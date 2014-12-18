//____________________________________________________________________________
/*!

\namespace genie::controls

\brief     Misc GENIE control constants

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           STFC, Rutherford Appleton Laboratory

\created   May 03, 2004

\cpright   Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _CONTROLS_H_
#define _CONTROLS_H_

namespace genie {
namespace controls {

// Maximum allowed number of iterations in rejection MC method
// before selecting a valid number
static const unsigned int kRjMaxIterations = 1000;

// Maximum allowed depth when GENIE is running in recursive mode
static const unsigned int kRecursiveModeMaxDepth = 100;

// Maximum allowed number of EVGThreadExceptions that is allowed
// to be caught by EventGenerator at a single event generation thread
static const unsigned int kMaxEVGThreadExceptions = 350;

// Default random number generator seed number. It can be overriden
// setting the $GSEED env. var. or by using RandomGen::SetSeed(int)
static const unsigned int kDefaultRandSeed = 65539;

static const double kASmallNum      = 1E-6;  
static const double kMinQ2Limit     = 1E-4;  // GeV^2
static const double kMinQ2Limit_VLE = 1E-10; // GeV^2
static const double kMinX           = 1E-4;
static const double kMinX_VHE       = 1E-5;
static const double kMaxX           = 1.-kASmallNum;
static const double kMinY           = 1E-4;
static const double kMinY_VHE       = 1E-5;
static const double kMaxY           = 1.-kASmallNum;

// KNO Hadronization model control parameters

// Default 'maximum' multiplicity for multiplicity probability distributions.
// This is not a 'hard limit'. If it is needed it will be extended internally
// by the KNO hadronization model.
static const int kMaxMultiplicity = 35;  

// Maximum number of attempts by the KNO hadronizer for finding a valid f/s
// hadronic system before going in error and quiting
static const unsigned int kMaxKNOHadSystIterations = 400;  

// Maximum number of attempts before producing an unweighted decay using the
// TGenPhaseSpace phase space generator
static const unsigned int kMaxUnweightDecayIterations = 1000;  

// Ma-like parameter used in variable transformations taking out the dipole 
// form factor form speeding up kinematical selection for QEL and RES events
static const double kMQD2 = 0.7;  

} // namespace controls
} // namespace genie

#endif // _CONTROLS_H_


