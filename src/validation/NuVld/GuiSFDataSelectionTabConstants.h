//_____________________________________________________________________________
/*!

\namespace  genie::nuvld::constants

\brief      NuValidator constants

\author     Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created    January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _SF_DATA_SELECTION_TAB_CONSTANTS_H_
#define _SF_DATA_SELECTION_TAB_CONSTANTS_H_

namespace genie {
namespace nuvld {
namespace constants {

// Constants used in the Structure Functions Selection Tab

static const char * kSFErrType[] = { "no error", "stat. only", "stat.+syst.", 0 };

static const char * kSFExperimentName[] =
{
   "NUTEV", "BCDMS", "CCFR", "CDHS", "EMC", "NMC", "SLAC", "WA59", 0
};
static const char * kSFTarget[] =
{
   "Proton", "Iron", "Neon", 0
};
static const char * kSFProbe[] =
{
   "nu_mu", "muon", "electron", 0
};
static const char * kSFName[] =
{
   "F2", "xF3", 0
};
static const char * kSFPlotVar[] =
{
   "Q2", "x", 0
};
static const char * kSFR[] =
{
   "QCD", "Whitlow", "Unknown", 0
};

} // namespace constants
} // namespace nuvld
} // namespace genie

#endif


