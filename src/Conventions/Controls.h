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
static const unsigned int kRecursiveModeMaxDepth = 20;

// maximum allowed number of EVGThreadExceptions that is allowed
// to be caught by EventGenerator at a single event generation thread
static const unsigned int kMaxEVGThreadExceptions = 5;

// Default random number generator seed number. It can be overriden
// setting the $GSEED env. var. or by using RandomGen::SetSeed(int)
static const unsigned int kDefaultRandSeed = 65539;

// Names of environmental variables to keep track of when saving
// the job configuration in TFolders along with the output event tree
static const char * kMCEnv[] =
{
  "GENIE", 
  "GEVGL", 
  "GSPLOAD", 
  "GSPSAVE", 
  "GMSGCONF",
  "ROOTSYS", 
  "CERNLIB", 
  "PYTHIA6", 
  "NEUGEN3PATH",
  "LIBXML2_INC", 
  "LIBXML2_LIB", 
  "LOG4CPP_INC", 
  "LOG4CPP_LIB",
  "CLHEP_INC", 
  "CLHEP_LIB", 
  "GEANT4_INC", 
  "GEANT4_LIB",
  "GPROFILER_LIB", 
  "DOXYGEN",
  "GOPT_ENABLE_NUVALIDATOR", 
  "GOPT_ENABLE_NEUGEN",
  "GOPT_ENABLE_GEANT_INTERFACE",
  "GOPT_ENABLE_FLUX_DRIVERS", 
  "GOPT_ENABLE_GEOM_DRIVERS",
  "GOPT_ENABLE_PROFILER", 
  "GOPT_ENABLE_DOXYGEN_DOC"
  "PATH", 
  "LD_LIBRARY_PATH", 
  0
};


} // namespace controls
} // namespace genie

#endif // _CONTROLS_H_


