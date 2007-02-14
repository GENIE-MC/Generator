
//____________________________________________________________________________
/*!

Names of environmental variables to keep track of when saving the job config 
in TFolders along with the output event tree

\author    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
           CCLRC, Rutherford Appleton Laboratory

\created   May 03, 2004

\cpright   Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
           All rights reserved.
           For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _ENV_SNAPSHOT_H_
#define _ENV_SNAPSHOT_H_

namespace genie {
namespace controls {

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
  "GOPT_ENABLE_DOXYGEN_DOC",
  "PATH", 
  "LD_LIBRARY_PATH", 
  0
};


} // namespace controls
} // namespace genie

#endif // _ENV_SNAPSHOT_H_


