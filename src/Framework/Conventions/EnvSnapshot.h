
//____________________________________________________________________________
/*!

Names of environmental variables to keep track of when saving the job config
in TFolders along with the output event tree

\author    Costas Andreopoulos <c.andreopoulos \at cern.ch>
           University of Liverpool

\created   May 03, 2004

\cpright   Copyright (c) 2003-2025, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org         
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
  "GSEED",
  "GSPLOAD",
  "GSPSAVE",
  "GMSGCONF",
  "GPRODMODE",
  "GALGCONF",
  "GCACHEFILE",
  "GUSERPHYSOPT",
  "GUNPHYSMASK",
  "ROOTSYS",
  "CERNLIB",
  "LHAPATH",
  "PYTHIA6",
  "LIBXML2_INC",
  "LIBXML2_LIB",
  "LOG4CPP_INC",
  "LOG4CPP_LIB",
  "GOPT_ENABLE_FLUX_DRIVERS",
  "GOPT_ENABLE_GEOM_DRIVERS",
  "PATH",
  "LD_LIBRARY_PATH",
  0
};

} // namespace controls
} // namespace genie

#endif // _ENV_SNAPSHOT_H_
