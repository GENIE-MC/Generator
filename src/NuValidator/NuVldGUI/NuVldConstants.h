//_____________________________________________________________________________
/*!

\namespace  genie::nuvld::constants

\brief      NuValidator constants

\author     Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created    January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NUVLD_CONSTANTS_H_
#define _NUVLD_CONSTANTS_H_

namespace genie {
namespace nuvld {
namespace constants {

//static const char * kSqlFileExtensions[] =
//{
// "All files", "*", "SQL files", "*.sql", 0, 0
//};

//static const char * kRootFileExtensions[] =
//{
//  "All files", "*", "ROOT files", "*.root", 0, 0
//};

static const char * kSavedPlotExtensions[] =
{
  "GIF", "*.gif", "ROOT macro", "*.C", "Encapsulated postscript", "*.eps", 0, 0
};

static const char * kFitters[] =
{
   "NONE", "SIMPLE", "SIMPLE/UWF", "NORM-FLOAT", "NORM-FLOAT/UWF", 0, 
};

static const char * kExperimentName[] =
{
   "Gargamelle", "BEBC", "ANL-12ft", "BNL-7ft", "FNAL-15ft", "SKAT",
   "CDHS", "CCFR", "CCFRR", "LSND", "IHEP-ITEP", "IHEP-JINR", "SERP-A1", "CHARM", 0
};

static const char * kExperimentMySQLName[] =
{
   "Gargamelle", "BEBC", "ANL_12FT", "BNL_7FT", "FNAL_15FT", "SKAT",
   "CDHS", "CCFR", "CCFRR", "LSND", "IHEP_ITEP", "IHEP_JINR", "SERP_A1", "CHARM", 0
};

static const char * kProcName[] =
{
   "Quasi Elastic", "Total", "Single Pion", "Multi Pion", "Coherent", 0
};

static const char * kProcMySQLName[] =
{
   "QES_XSEC", "TOT_XSEC", "SPP_XSEC", "MPP_XSEC", "COH_XSEC", 0
};

static const char * kNuType[] =
{
   "nu_e", "nu_e_bar", "nu_mu", "nu_mu_bar", "nu_tau", "nu_tau_bar", 0
};

static const char * kNuTypeMySQLName[] =
{
   "nu_e", "nu_e_bar", "nu_mu", "nu_mu_bar", "nu_tau", "nu_tau_bar", 0
};

static const char * kTarget[] =
{
   "Hydrogen", "Deuterium", "Neon", "Iron", "Aluminium", "Carbon", "Propane", "Freon", 0
};

static const char * kTargetMySQLName[] =
{
   "Hydrogen", "Deuterium", "Neon", "Iron", "Aluminium", "Carbon", "Propane", "Freon", 0
};

static const char * kXSecErrType[] = { "no error", "stat. only", "stat.+syst.", 0 };

static const char * kXSecErrDrawOpt[] = { "noXsec", "stat", "all", 0 };

static const double kEmin =   0.1;
static const double kEmax = 120.00;

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


