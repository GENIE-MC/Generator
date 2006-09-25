//____________________________________________________________________________
/*!

\namespace genie::constants::bba2005

\brief     Coefficients for the BBA2003 (Bradford,Budd,Bodek,Arrington) elastic 
           form factor parameterization.
           Note: these values will be used as 'defaults' if no other values
           are given at the BBA05 model XML config file.

\ref       R.Bradford, NuINT05 proceedings, hep-ex/0602017

\author    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
           CCLRC, Rutherford Appleton Laboratory

\created   September 25, 2006

\cpright   Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
           All rights reserved.
           For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _BBA2005_CONSTANTS_H_
#define _BBA2005_CONSTANTS_H_

namespace genie {
namespace constants {
namespace bba2005 {

static const double kGep_a0 =    1.;
static const double kGep_a1 =   -0.0578; // +/- 0.1660
static const double kGep_a2 =    0.;
static const double kGep_b1 =   11.100;  // +/- 0.217
static const double kGep_b2 =   13.60;   // +/- 1.39
static const double kGep_b3 =   33.00;   // +/- 8.95
static const double kGep_b4 =    0.;

static const double kGmp_a0 =    1.;
static const double kGmp_a1 =    0.1500; // +/- 0.0312
static const double kGmp_a2 =    0.;
static const double kGmp_b1 =   11.100;  // +/- 0.103
static const double kGmp_b2 =   19.600;  // +/- 0.281
static const double kGmp_b3 =    7.540;  // +/- 0.967
static const double kGmp_b4 =    0.;

static const double kGen_a0 =    0.;
static const double kGen_a1 =    1.250; // +/-   0.368
static const double kGen_a2 =    1.30;  // +/-   1.99
static const double kGen_b1 =   -9.86;  // +/-   6.46
static const double kGen_b2 =  305.0;   // +/-  28.6
static const double kGen_b3 = -758.0;   // +/-  77.5
static const double kGen_b4 =  802.;    // +/- 156.

static const double kGmn_a0 =    1.;
static const double kGmn_a1 =    1.810; // +/- 0.402
static const double kGmn_a2 =    0.;
static const double kGmn_b1 =   14.100; // +/- 0.597
static const double kGmn_b2 =   20.70;  // +/- 2.55
static const double kGmn_b3 =   68.7;   // +/- 14.1
static const double kGmn_b4 =    0;

// Gep/Gmp assummed const for Q2 > Q2Max
//static const double kQ2Max   =  6.0;

} // namespace bba2005
} // namespace constants
} // namespace genie

#endif // _BBA2005_CONSTANTS_H_



