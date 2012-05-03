//_____________________________________________________________________________
/*!

\namespace  genie::nuvld::constants

\brief      NuValidator constants

\author     Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created    January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _ELECTRON_DATA_SELECTION_TAB_CONSTANTS_H_
#define _ELECTRON_DATA_SELECTION_TAB_CONSTANTS_H_

namespace genie {
namespace nuvld {
namespace constants {

// Constants used in the eN data Selection Tab

static const char * kElExperiment[] =
{
   "E133", "E140", "E140X", "E49A10", "E49A6", "E49B", "E61", "E87", 
   "E891", "E8920", "JLAB", "NE11", "ONEN1HAF", 0
};

static const char * kElTarget[] =
{
   "Hydrogen", "Deuterium", 0
};

static const int kNElVarRangeFrames = 9;

static const char * kElVarFrameName[kNElVarRangeFrames] =
{
  "E (GeV)", "Ep (GeV)", "Theta (degrees)", "Q^2 (GeV^2)",
  "W^2 (GeV^2)", "v (GeV)", "x", "Epsilon", "Gamma"
};

static const char * kElVarMySQLName[kNElVarRangeFrames] =
{
  "E", "EP", "Theta", "Q2", "W2", "Nu", "x", "Epsilon", "Gamma"
};

static const double kElVarMin[kNElVarRangeFrames] =
{
  0, 0, 0, 0, 0, 0, 0, 0, 0
};

static double kElVarMax[kNElVarRangeFrames] =
{
  100, 100, 100, 100, 100, 100, 100, 100, 100
};

} // namespace constants
} // namespace nuvld
} // namespace genie

#endif


