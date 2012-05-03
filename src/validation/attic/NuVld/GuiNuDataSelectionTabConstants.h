//_____________________________________________________________________________
/*!

\namespace  genie::nuvld::constants

\brief      NuValidator constants

\author     Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created    January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NUVLD_CONSTANTS_1_H_
#define _NUVLD_CONSTANTS_1_H_

namespace genie {
namespace nuvld {
namespace constants {

// Constants used in the vN data Selection Tab

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

static const char * kXSecErrType[] = { "none", "stat", "syst", "stat+syst", 0 };

static const double kEmin =   0.1;
static const double kEmax = 120.00;

} // namespace constants
} // namespace nuvld
} // namespace genie

#endif


